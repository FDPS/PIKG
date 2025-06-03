# MN-Core2 params
num_group = 4
num_l2b_per_group = 2
num_l1b_per_l2b = 8
num_mab_per_l1b = 16
num_pe_per_mab = 4

num_l2b = num_group * num_l2b_per_group
num_l1b = num_l2b * num_l1b_per_l2b
num_mab = num_l1b * num_mab_per_l1b
num_pe = num_mab * num_pe_per_mab

num_pe_per_l1b = num_mab_per_l1b * num_pe_per_mab
num_pe_per_l2b = num_l1b_per_l2b * num_pe_per_l1b

cycle_per_step = 4

l2bm_chunk_size = 64

HW_BIT_SIZE = 16
SW_BIT_SIZE = 32
LW_BIT_SIZE = 64
SW_IN_LW = LW_BIT_SIZE//SW_BIT_SIZE
HW_IN_LW = LW_BIT_SIZE//HW_BIT_SIZE

LMEM_SIZE = 4096 # in short word
GRF_SIZE = 512 # in short word
DRAM_SIZE = 4*1024*1024*1024 // 8 # in long word

LMEM_BANK_OFFSET=256
GRF_BANK_OFFSET=256
DRAM_OFFSET_SIZE=4*1024*1024//8
DRAM_IN_OFFSET = 0
DRAM_OUT_OFFSET = DRAM_IN_OFFSET + DRAM_OFFSET_SIZE
DRAM_IN_J_OFFSET = DRAM_OUT_OFFSET + DRAM_OFFSET_SIZE
DRAM_KERNEL_OFFSET = DRAM_IN_J_OFFSET + 2*DRAM_OFFSET_SIZE

MAX_KERNEL_LAUNCH = 16
