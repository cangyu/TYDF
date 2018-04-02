#include <iostream>
#include <vector>

using namespace std;

typedef enum {WALL, ONE_TO_ONE, FAR} NMF_BCType;

class NMF_BCRange
{
private:
	uint32_t blk_seq;
	uint32_t face_seq;
	uint32_t pri_start_seq, pri_end_seq;
	uint32_t sec_start_seq, sec_end_seq;

public:
	NMF_BCRange():blk_seq(0), face_seq(0), pri_start_seq(0), pri_end_seq(0), sec_start_seq(0), sec_end_seq(0) {}
	NMF_BCRange(uint32_t b, uint32_t f, uint32_t s1, uint32_t e1, uint32_t s2, uint32_t e2) :blk_seq(b), face_seq(f), pri_start_seq(s1), pri_end_seq(e1), sec_start_seq(s2), sec_end_seq(e2) {}
};

class NMF_Entry
{
private:
	NMF_BCType bc;
	NMF_BCRange *rg1, *rg2;
	bool swap;

public:
	NMF_Entry(NMF_BCType t):bc(t)
	{
		rg1 = nullptr;
		rg2 = nullptr;
		swap = false;
	}

} ;

typedef struct
{
	uint32_t i_dim, j_dim, k_dim;
}NMF_BLKDim;

typedef struct
{
	uint32_t num;
	NMF_BLKDim *dim;
	NMF_Entry *entry;
} NMF;

void read_nmf(FILE *fin, NMF *dst)
{
	for (int i = 0; i < 4; i++)
		fscanf(fin, "%*[^\n]%*c");

	fscanf(fin, "%ud", &dst->num);
}

int main(int argc, char **argv)
{
	if (argc != 3)
	{
		printf("Usage: gridglue grid.fmt map.nmf\n");
		return -1;
	}

	//FILE *gfp = fopen(argv[1], "r");
	FILE *mfp = fopen(argv[2], "r");

	NMF mapping;
	read_nmf(mfp, &mapping);

	printf("%d\n", mapping.num);


	//fclose(gfp);
	fclose(mfp);

	return 0;
}