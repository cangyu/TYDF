'''Grid parameters'''
la = self.profile[0].chord_len
lt = 30 * la
r = 10 * la
inner_spn = self.profile[-1].z_offset
outer_spn = 20 * inner_spn








blk_param_list = []

b0_tfi_grid = 
blk_list.append(b0_tfi_grid)
blk_param_list.append([knot_dist[3], knot_dist[0], knot_dist[2]])

b1_tfi_grid = 
blk_list.append(b1_tfi_grid)
blk_param_list.append([knot_dist[3], knot_dist[7], knot_dist[2]])

b2_tfi_grid = 
blk_list.append(b2_tfi_grid)
blk_param_list.append([knot_dist[0], knot_dist[3], knot_dist[2]])

b3_tfi_grid = 
blk_list.append(b3_tfi_grid)
blk_param_list.append([knot_dist[3], knot_dist[0], knot_dist[4]])

b4_tfi_grid = 
blk_list.append(b4_tfi_grid)
blk_param_list.append([knot_dist[3], knot_dist[7], knot_dist[4]])

b5_tfi_grid = 
blk_list.append(b5_tfi_grid)
blk_param_list.append([knot_dist[0], knot_dist[3], knot_dist[4]])

blk_param_list.append([knot_dist[0], knot_dist[1], knot_dist[2]])
blk_param_list.append([knot_dist[5], knot_dist[0], knot_dist[2]])
blk_param_list.append([knot_dist[6], knot_dist[0], knot_dist[2]])
blk_param_list.append([knot_dist[0], knot_dist[1], knot_dist[4]])
blk_param_list.append([knot_dist[5], knot_dist[0], knot_dist[4]])
blk_param_list.append([knot_dist[6], knot_dist[0], knot_dist[4]])
blk_param_list.append([knot_dist[7], knot_dist[6], knot_dist[4]])

def report(msg):
    print('Process {} : {}'.format(os.getpid(), msg))



'''����, �߽�����, �ڽӹ�ϵ'''
blk = [b0_tfi_grid.grid,
        b1_tfi_grid.grid,
        b2_tfi_grid.grid,
        b3_tfi_grid.grid,
        b4_tfi_grid.grid,
        b5_tfi_grid.grid,
        b6_tfi_grid.grid,
        b7_tfi_grid.grid,
        b8_tfi_grid.grid,
        b9_tfi_grid.grid,
        b10_tfi_grid.grid,
        b11_tfi_grid.grid,
        b12_tfi_grid.grid]

bc = [(BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b0
        (BCType.Wall, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Symmetry, BCType.Interior),  # b1
        (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b2
        (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b3
        (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField),  # b4
        (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b5
        (BCType.Wall, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Symmetry, BCType.Interior),  # b6
        (BCType.Interior, BCType.Interior, BCType.Wall, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b7
        (BCType.Interior, BCType.Interior, BCType.Wall, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b8
        (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField),  # b9
        (BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b10
        (BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b11
        (BCType.Interior, BCType.Interior, BCType.Interior, BCType.Interior, BCType.Wall, BCType.PressureFarField)  # b12
        ]

adj = [((6, 3), (0, 1), 1, True),
        ((0, 2), (0, 0), 0, False),
        ((1, 4), (0, 3), 1, False),
        ((0, 4), (0, 0), 0, False),
        ((0, 0), (0, 5), 1, False),
        ((0, 6), (3, 5), 0, False),

        ((0, 0), (1, 1), 1, False),
        ((1, 2), (0, 0), 0, False),
        ((2, 1), (1, 3), 1, True),
        ((0, 0), (1, 5), 1, False),
        ((1, 6), (4, 5), 0, False),

        ((2, 2), (0, 0), 0, False),
        ((8, 1), (2, 3), 1, True),
        ((2, 4), (0, 0), 0, False),
        ((0, 0), (2, 5), 1, False),
        ((2, 6), (5, 5), 0, False),

        ((9, 3), (3, 1), 1, True),
        ((3, 2), (0, 0), 0, False),
        ((4, 4), (3, 3), 1, False),
        ((3, 4), (0, 0), 0, False),
        ((3, 6), (0, 0), 0, False),

        ((12, 3), (4, 1), 1, True),
        ((4, 2), (0, 0), 0, False),
        ((5, 1), (4, 3), 1, True),
        ((4, 6), (0, 0), 0, False),

        ((5, 2), (0, 0), 0, False),
        ((11, 1), (5, 3), 1, True),
        ((5, 4), (0, 0), 0, False),
        ((5, 6), (0, 0), 0, False),

        ((0, 0), (6, 1), 1, False),
        ((6, 2), (0, 0), 0, False),
        ((6, 4), (7, 2), 0, True),
        ((0, 0), (6, 5), 1, False),
        ((6, 6), (9, 5), 0, False),

        ((8, 2), (7, 1), 1, False),
        ((0, 0), (7, 3), 1, False),
        ((7, 4), (0, 0), 0, False),
        ((0, 0), (7, 5), 1, False),
        ((7, 6), (10, 5), 0, False),

        ((0, 0), (8, 3), 1, False),
        ((8, 4), (0, 0), 0, False),
        ((0, 0), (8, 5), 1, False),
        ((8, 6), (11, 5), 0, False),

        ((12, 2), (9, 1), 1, False),
        ((9, 2), (0, 0), 0, False),
        ((9, 4), (10, 2), 0, True),
        ((9, 6), (0, 0), 0, False),

        ((11, 2), (10, 1), 1, False),
        ((12, 4), (10, 3), 1, False),
        ((10, 4), (0, 0), 0, False),
        ((10, 6), (0, 0), 0, False),

        ((12, 1), (11, 3), 1, True),
        ((11, 4), (0, 0), 0, False),
        ((11, 6), (0, 0), 0, False),

        ((0, 0), (12, 5), 1, False),
        ((12, 6), (0, 0), 0, False)]

'''����MSH�ļ�'''
msh = XF_MSH.from_str3d_multi(blk, bc, adj)
