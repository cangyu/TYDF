#include "nmf.h"

int main(int argc, char **argv)
{
	/*
	if (argc != 3)
	{
		cout << "Usage: gridglue map.nmf grid.fmt" << endl;
		return -1;
	}
	*/

	NMF bc_mp("../../test/mapping2.nmf");

	uint32_t total_cell_num = 0, total_face_num = 0, total_node_num = 0;
	for (uint32_t i = 0; i < bc_mp.blk_num; i++)
	{
		total_cell_num += bc_mp.blk_dim[i].cell_num();
		total_face_num += bc_mp.blk_dim[i].face_num();
		total_node_num += bc_mp.blk_dim[i].internal_node_num();
	}

	for (uint32_t i = 0; i < bc_mp.mapping_entry.size(); i++)
		if (bc_mp.mapping_entry[i]->bc == bc_o2o)
			total_face_num -= bc_mp.mapping_entry[i]->rg1->face_num();

	vector<vector<uint32_t>> shell_node_seq(bc_mp.blk_num);
	for (uint32_t i = 0; i < bc_mp.blk_num; i++)
	{
		shell_node_seq[i].resize(bc_mp.blk_dim[i].shell_node_num());
		memset(shell_node_seq[i].data(), 0, shell_node_seq[i].size() * sizeof(uint32_t));
	}

	uint32_t cnt = total_node_num + 1;
	for (uint32_t k = 0; k < bc_mp.mapping_entry.size(); k++)
	{
		if (bc_mp.mapping_entry[k]->bc == bc_o2o)
		{
			uint32_t blk = bc_mp.mapping_entry[k]->rg1->blk_seq - 1;
			uint32_t face = bc_mp.mapping_entry[k]->rg1->face_seq;

			for (uint32_t pri = bc_mp.mapping_entry[k]->rg1->pri_start_seq; pri <= bc_mp.mapping_entry[k]->rg1->pri_end_seq; pri++)
			{
				for (uint32_t sec = bc_mp.mapping_entry[k]->rg1->sec_start_seq; sec <= bc_mp.mapping_entry[k]->rg1->sec_end_seq; sec++)
				{
					uint32_t ii = 0, jj = 0, kk = 0;
					bc_mp.blk_dim[blk].node_coordinate(face, pri, sec, ii, jj, kk);
					uint32_t s_idx = bc_mp.blk_dim[blk].shell_node_idx(ii, jj, kk);
					if (shell_node_seq[blk][s_idx] == 0)
					{
						vector<pair<uint32_t, uint32_t>> g;
						g.push_back(make_pair(blk, s_idx));
						bc_mp.find_all_occurrence(blk, ii, jj, kk, g);
						for (uint32_t t = 0; t < g.size(); t++)
						{
							uint32_t b = g[t].first;
							uint32_t idx = g[t].second;
							shell_node_seq[b][idx] = cnt;
						}
						++cnt;
					}
				}
			}
		}
	}

	for (uint32_t k = 0; k < bc_mp.mapping_entry.size(); k++)
	{
		if (bc_mp.mapping_entry[k]->bc != bc_o2o)
		{
			uint32_t blk = bc_mp.mapping_entry[k]->rg1->blk_seq - 1;
			uint32_t face = bc_mp.mapping_entry[k]->rg1->face_seq;

			for (uint32_t pri = bc_mp.mapping_entry[k]->rg1->pri_start_seq; pri <= bc_mp.mapping_entry[k]->rg1->pri_end_seq; pri++)
			{
				for (uint32_t sec = bc_mp.mapping_entry[k]->rg1->sec_start_seq; sec <= bc_mp.mapping_entry[k]->rg1->sec_end_seq; sec++)
				{
					uint32_t ii = 0, jj = 0, kk = 0;
					bc_mp.blk_dim[blk].node_coordinate(face, pri, sec, ii, jj, kk);
					uint32_t s_idx = bc_mp.blk_dim[blk].shell_node_idx(ii, jj, kk);
					if (shell_node_seq[blk][s_idx] == 0)
						shell_node_seq[blk][s_idx] = cnt++;
				}
			}
		}
	}

	total_node_num = cnt - 1;

	cout << "Cell num: " << total_cell_num << endl;
	cout << "Face num: " << total_face_num << endl;
	cout << "Node num: " << total_node_num << endl;

	bc_mp.save("../../test/test2.nmf");

	return 0;
}