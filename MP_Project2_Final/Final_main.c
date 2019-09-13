#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#define SIZE 262144
#define MAX_LINES 131072  // Total cache of size  in our project 256KB
#define MAX_WAYS 16

uint32_t LRU[MAX_LINES][MAX_WAYS];
bool d[MAX_LINES][MAX_WAYS];
bool valid[MAX_LINES][MAX_WAYS];
uint32_t tag[MAX_LINES][MAX_WAYS];
uint8_t  N, bl, w = 24;
uint32_t read_mem_to_cache_count;
uint32_t write_mem_read_count;
uint32_t read_hit_count, write_cache_miss_count, read_miss_count, read_replace_count, write_replace_count, read_memory_count, write_memory_count;
uint32_t write_mem_dirty_count, read_writeback_count, write_writeback_count, write_cache_count, read_cache_count, write_through_count, read_cache_to_CPU_count = 0;
uint8_t WS; // write strategy -> WS=1--> writeback, WS=2--> WTA, WS=3--> WTNA
uint32_t s = 262144;
uint32_t t, l, b, L, B;
uint32_t line, tag_current;
uint32_t temp_LRU, LRU_Max, way1, line1, total_bytes_read;
uint32_t write_cache_hit_count;
uint32_t total_bytes_write;
int32_t way;
bool hit;
uint32_t n = 256;
bool hit;
static double p[256];
uint32_t row = 256;
uint32_t column = 256;
static double a[256][256];
static double var[256][256];
uint8_t zerocache = 0;
static double x[256];
uint32_t flush_count;

void zeroCache()
{
	for (line = 0; line < MAX_LINES; line++)
	{
		for (way = 0; way < MAX_WAYS; way++)
		{
			valid[line][way] = 0;
			d[line][way] = 0;
			LRU[line][way] = way;

		}
	}
	zerocache++;
	printf("ZeroCache: %d\n", zerocache);
}

void flush_cache()
{
	uint32_t i;
	uint8_t j;
	for (i = 0; i < L; i++)
	{
		for (j = 0; j < N; j++)
		{
			if (d[i][j] == 1 && WS == 1)
			{
				flush_count++;
				d[i][j] = 0;
			}
		}
	}
}
void clear_performance_counters()
{
	write_cache_count = 0;
	write_mem_dirty_count = 0;
	write_mem_read_count = 0;
	write_cache_miss_count = 0;
	read_cache_count = 0;
	read_writeback_count = 0;
	read_hit_count = 0;
	read_mem_to_cache_count = 0;
	read_cache_to_CPU_count = 0;
	read_miss_count = 0;
	write_through_count = 0;
	write_cache_hit_count = 0;
	read_replace_count = 0;
	write_replace_count = 0;
	write_memory_count = 0;
	read_memory_count = 0;
	total_bytes_write = 0;
	total_bytes_read = 0;
	write_writeback_count = 0;
	flush_count = 0;

}

void calc(uint8_t N, uint8_t bl)
{
	B = bl * 2;
	L = s / (N* B);
	b = (abs)(log(B) / log(2));
	l = (abs)(log(L) / log(2));
}

uint32_t getLine(uint32_t add, uint8_t bl, uint8_t N)
{
	uint32_t i, j;
	i = (add >> b);
	j = (pow(2, l) - 1);
	line = i & j;
	return (line);
}
uint32_t getTag(uint32_t add, uint8_t bl, uint8_t N)
{
	uint32_t t;
	tag_current = (add >> (b + l));
	t = w - b - l;
	tag_current = tag_current & t;
	return tag_current;
}
bool isdirty(uint32_t line, uint32_t way)
{
	if (d[line][way] == 1 && valid[line][way] == 1)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void invalidate(uint32_t line, uint32_t way)
{
	valid[line][way] = 0;
}
void validate(uint32_t line, uint32_t way)
{
	valid[line][way] = 1;
}
uint32_t find_oldestLRU(uint32_t line, uint8_t N)
{
	uint32_t i, LRUway = 0, tempLRU;
	for (i = 0; i < N; i++)
	{
		if (LRU[line][i] == N - 1)
		{
			tempLRU = LRU[line][i];
			LRUway = i;
		}
	}
	LRU_Max = temp_LRU;
	return LRUway;
}

void update_LRU(uint32_t line, uint32_t way)
{
	uint32_t i;
	for (i = 0; i < N; i++)
	{
		LRU[line][i] < LRU_Max;
		LRU[line][i]++;
	}
	LRU[line][way] = 0;
}
void readcache_to_CPU()
{
	read_cache_to_CPU_count++;
}
void readmem_to_Cache()
{
	read_mem_to_cache_count++;
}
void set_Tag(uint32_t line, uint32_t way, uint32_t tag_current)
{
	tag[line][way] = tag_current;
}
int8_t findTag(uint32_t line, uint32_t tag_current, uint8_t N)
{
	uint32_t j;
	for (j = 0; j < N; j++)
	{
		if ((tag[line][j] == tag_current) && (valid[line][j] == 1))
		{
			way1 = j;
			return way1;
		}
	}
	return -1;
}

void read_writeback()
{
	read_writeback_count++;
}

void RC(void *add, uint8_t N, uint8_t bl)
{
	uint32_t addnew1 = (uint32_t)add;
	read_cache_count++;
	line = getLine(addnew1, bl, N);
	tag_current = getTag(addnew1, bl, N);
	way = findTag(line, tag_current, N);
	hit = (way != -1);
	if (!hit)
	{
		read_miss_count++;
		way = find_oldestLRU(line, N);
		if (valid[line][way] == 1)
		{
			read_replace_count++;
			if (isdirty(line, way))
			{
				read_writeback();
				d[line][way] = 0;
			}
		}
		invalidate(line, way); // make valid bits 0
		set_Tag(line, way, tag_current);
		line1 = line;
		way1 = way;
		update_LRU(line1, way1);
		readmem_to_Cache(); // increment read_cache counter
		validate(line, way);// make valid bits 1
	}
	else
	{
		read_hit_count++;
		update_LRU(line, way);
		readcache_to_CPU();
	}

}

void RM(void *pmem, uint32_t size)
{
	read_memory_count++;
	int32_t last_line = -1;
	uint32_t add = (uint32_t)pmem;
	for (uint32_t i = 0; i < size; i++)
	{
		line = getLine(add, bl, N);
		if (line != last_line)
		{
			RC(add, N, bl);
			last_line = line;
		}
		add++;
	}
	total_bytes_read += size;
}

void write_cache_block()
{

}

void write_mem_to_cache()
{
	write_mem_read_count++;
}

void WC(void *add, uint8_t N, uint8_t bl)
{
	uint32_t addnew = (uint32_t)add;
	write_cache_count++;
	tag_current = getTag(addnew, bl, N);
	line = getLine(addnew, bl, N);
	way = findTag(line, tag_current, N);
	hit = (way != -1);
	if (!hit)
	{
		write_cache_miss_count++;
	}
	if ((!hit) && (WS != 3))
	{
		way = find_oldestLRU(line, N);
		if (valid[line][way] == 1)
		{
			write_replace_count++;
			if (isdirty(line, way))
			{
				write_mem_dirty_count++;
				d[line][way] = 0;
			}
		}
		invalidate(line, way);
		write_cache_block();
		validate(line, way);
		set_Tag(line, way, tag_current);
		line1 = line;
		way1 = way;
	}

	if ((hit) || (WS != 3))
	{
		if (hit)
		{
			write_cache_hit_count++;
		}
		update_LRU(line1, way1);
		write_mem_to_cache();
	}


	if ((WS == 2) || (WS == 3))
	{
		write_through_count++;
	}
	if (WS == 1)
	{
		write_writeback_count++;
		d[line][way] = 1;
	}
}

void WM(void *pmem, uint32_t size)
{
	write_memory_count++;
	int32_t last_line = -1;
	uint32_t add = (uint32_t)pmem;
	for (uint32_t i = 0; i < size; i++)
	{
		line = getLine(add, bl, N);
		if (line != last_line)
		{
			WC(add, N, bl);
			last_line = line;
		}
		add++;
	}
	total_bytes_write += size;

}




void choldc(double a[256][256], uint32_t n, double p[])
{
	uint32_t i, j;
	int32_t k;
	double sum;

	WM(&i, sizeof(uint32_t));
	RM(&i, sizeof(uint32_t));
	RM(&n, sizeof(uint32_t));
	for (i = 0; i < n; i++)
	{
		WM(&j, sizeof(uint32_t));
		RM(&i, sizeof(uint32_t));
		RM(&j, sizeof(uint32_t));
		RM(&n, sizeof(uint32_t));
		for (j = i; j < n; j++)
		{
			RM(&i, sizeof(uint32_t));
			RM(&j, sizeof(uint32_t));
			RM(&a[i][j], sizeof(double));
			WM(&sum, sizeof(double));
			sum = a[i][j];

			WM(&k, sizeof(uint32_t));
			RM(&k, sizeof(uint32_t));
			RM(&i, sizeof(uint32_t));

			for (k = i - 1; k >= 0; k--)
			{
				RM(&i, sizeof(uint32_t));
				RM(&j, sizeof(uint32_t));
				RM(&k, sizeof(uint32_t));
				RM(&a[i][k], sizeof(double));
				RM(&a[j][k], sizeof(double));
				WM(&sum, sizeof(double));
				sum -= a[i][k] * a[j][k];

			}
			RM(&j, sizeof(uint32_t));
			RM(&i, sizeof(uint32_t));

			if (i == j) {
				RM(&sum, sizeof(double));
				if (sum <= 0.0)
				{
					printf("choldc failed");
					//	exit(1);
				}
				RM(&sum, sizeof(double));
				RM(&i, sizeof(uint32_t));
				WM(&p[i], sizeof(double));
				p[i] = sqrt(sum);

			}
			else
				RM(&sum, sizeof(double));
			RM(&i, sizeof(uint32_t));
			RM(&j, sizeof(uint32_t));
			RM(&p[i], sizeof(double));
			WM(&a[j][i], sizeof(double));
			a[j][i] = sum / p[i];

			RM(&j, sizeof(uint32_t));
			RM(&n, sizeof(uint32_t));
			WM(&j, sizeof(uint32_t));

		}
		RM(&i, sizeof(uint32_t));
		RM(&n, sizeof(uint32_t));
		WM(&i, sizeof(uint32_t));
	}

}


void cholsl(double a[256][256], int n, double p[], float b[], double x[])
{
	uint32_t i;
	int32_t k;
	double sum;
	WM(&i, sizeof(uint32_t));
	RM(&i, sizeof(uint32_t));
	RM(&n, sizeof(uint32_t));
	for (i = 0; i <= n; i++) {
		RM(&i, sizeof(uint32_t));
		RM(&b[i], sizeof(double));
		WM(&sum, sizeof(double));

		sum = b[i];
		WM(&k, sizeof(uint32_t));
		RM(&k, sizeof(uint32_t));
		RM(&k, sizeof(uint32_t));
		for (k = i - 1; k >= 1; k--)
		{
			RM(&i, sizeof(uint32_t));
			RM(&k, sizeof(int32_t));
			RM(&a[i][k], sizeof(double));
			RM(&x[k], sizeof(double));
			WM(&sum, sizeof(double));

			sum -= a[i][k] * x[k];

			WM(&k, sizeof(uint32_t));
			RM(&k, sizeof(uint32_t));
			RM(&k, sizeof(uint32_t));

		}
		RM(&i, sizeof(uint32_t));
		RM(&p[i], sizeof(double));
		RM(&sum, sizeof(double));
		WM(&x[i], sizeof(double));
		x[i] = sum / p[i];

		WM(&i, sizeof(uint32_t));
		RM(&i, sizeof(uint32_t));
		RM(&n, sizeof(uint32_t));
	}
	WM(&i, sizeof(uint32_t));
	RM(&i, sizeof(uint32_t));
	RM(&n, sizeof(uint32_t));
	for (i = n; i >= 1; i--) {
		RM(&i, sizeof(uint32_t));
		RM(&x[i], sizeof(uint32_t));
		WM(&i, sizeof(uint32_t));
		sum = x[i];
		WM(&k, sizeof(uint32_t));
		RM(&k, sizeof(uint32_t));
		RM(&n, sizeof(uint32_t));
		for (k = i + 1; k <= n; k++)
		{

			RM(&i, sizeof(uint32_t));
			RM(&k, sizeof(int32_t));
			RM(&a[k][i], sizeof(double));
			RM(&x[k], sizeof(double));
			WM(&sum, sizeof(double));
			sum -= a[k][i] * x[k];

			WM(&k, sizeof(uint32_t));
			RM(&k, sizeof(uint32_t));
			RM(&n, sizeof(uint32_t));

		}
		RM(&i, sizeof(uint32_t));
		RM(&p[i], sizeof(double));
		RM(&sum, sizeof(double));
		WM(&x[i], sizeof(double));

		x[i] = sum / p[i];
		WM(&i, sizeof(uint32_t));
		RM(&i, sizeof(uint32_t));
		RM(&n, sizeof(uint32_t));
	}

}


int file_write(FILE *file)
{
#pragma warning(disable:4996)
	if (file == NULL)
	{
		return -1;
	}

	fprintf(file, "%d\t", WS);
	fprintf(file, "%d\t", N);
	fprintf(file, "%d\t", bl);
	fprintf(file, "%d\t", read_memory_count);
	fprintf(file, "%d\t", read_cache_count);
	fprintf(file, "%d\t", read_miss_count);
	fprintf(file, "%d\t", read_hit_count);
	fprintf(file, "%d\t", read_mem_to_cache_count);
	fprintf(file, "%d\t", read_writeback_count);
	fprintf(file, "%d\t", read_replace_count);
	fprintf(file, "%d\t", total_bytes_read);
	fprintf(file, "%d\t", write_memory_count);
	fprintf(file, "%d\t", write_cache_count);
	fprintf(file, "%d\t", write_cache_hit_count);
	fprintf(file, "%d\t", write_cache_miss_count);
	fprintf(file, "%d\t", write_writeback_count);
	fprintf(file, "%d\t", write_through_count);
	fprintf(file, "%d\t", write_replace_count);
	fprintf(file, "%d\t", write_mem_dirty_count);
	fprintf(file, "%d\t", total_bytes_write);
	fprintf(file, "%d\t", flush_count);
	fprintf(file, "\n");

	return 0;
}


int main() {
#pragma align 16;
	FILE *fp;
	fp = fopen("OutputFile.xls", "w+");
	fprintf(fp, "WS\tN\tbl\tReadMemoryCount\tReadCacheCount\tReadMissCount\tReadHitCount\tReadMemortToCacheCount\tReadWriteBackCount\tReadReplaceCount\tTotalBytesRead\tWriteMemoryCount\tWriteCacheCount\tWriteHitCount\tWriteMissCount\tWriteWriteBackCount\tWriteWriteThroughCount\tWriteReplaceCount\tWriteMemoryDirtyCount\tTotalBytesWritten\tFlushCount\n");
	double b[256];
	//test func
	//uint8_t mat[262144];
	//uint64_t sum = 0;
	for (WS = 1; WS <= 3; WS++)
	{
		for (N = 1; N <= 16; N *= 2)
		{
			for (bl = 1; bl <= 8; bl *= 2)
			{
				calc(N, bl);
				zeroCache();
				clear_performance_counters();

				//N = 4; bl = 2;
				for (int i = 0; i < 256; i++)
				{
					for (int j = 0; j < 256; j++)
					{
						a[i][j] = (rand() % 10);
						if (i == j)
						{
							var[i][j] = 1;
						}
						else
						{
							var[i][j] = 0;
						}
					}
				}
				for (int i = 0; i < 256; i++)
				{
					for (int j = 0; j < 256; j++)
					{
						if (j < i)
						{
							a[i][j] = 0.5 *(a[i][j] + a[j][i]);
							a[j][i] = a[i][j];
						}

					}
				}
				for (int i = 0; i < 256; i++)
				{
					for (int j = 0; j < 256; j++)
					{

						a[i][j] += (2 * n) * var[i][j];


					}

				}
				for (int j = 0; j < 256; j++)
				{
					b[j] = (rand() % 10);
				}

				choldc(a, n, p);
				cholsl(a, n, p, b, x);

				//test
				/*printf("Read mem to cache count is %d \n", read_mem_to_cache_count);
				printf("Read hit count is %d \n", read_hit_count);
				printf("Read miss countis %d\n", read_miss_count);
				printf("Read writeback count is %d\n", read_writeback_count);
				printf("total bytes read is %d\n", total_bytes_read);
				printf("read cache count %d\n", read_cache_count);
				printf("read replace count %d\n", read_replace_count);
				printf("read memory count %d\n", read_memory_count);
				printf("WS = %d\n", WS);
				printf("write mem dirty count %d\n", write_mem_dirty_count);
				printf("write cache miss count %d\n", write_cache_miss_count);
				printf("write cache hit count %d\n", write_cache_hit_count);
				printf("write through count %d\n", write_through_count);
				printf("total bytes write is %d\n", total_bytes_write);
				printf("write cache count %d\n", write_cache_count);
				printf("write replace count %d\n", write_replace_count);
				printf("write writeback count %d\n", write_writeback_count);
				printf("write memory count %d\n", write_memory_count);*/

				if (WS == 1)
				{
					flush_cache();
				}

				file_write(fp);

			}
		}
	}
	//test
	/*printf("Print the decomposed matrix \n");
	for (int k = 0; k < row; k++) {
		for (int j = 0; j < column; j++) {
			printf(" %f \t", a[k][j]);
		}
		printf("\n");
	}

	printf("\nDiagonal Element matrix\n");
	for (int k = 0; k < row; k++) {
		printf(" %f \t", p[k]);
	}*/

	fclose(fp);
	return 0;


}