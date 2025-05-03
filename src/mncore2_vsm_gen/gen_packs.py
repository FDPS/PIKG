import sys
import json
import struct
import copy

import os

from codehandler import *
from generate_core_instlist import get_reg_name, get_reg_id, get_reg_type

from constants import *

exe=sys.argv.pop(0)
max_subwalk_size = 0
while len(sys.argv) > 0:
    opt = sys.argv.pop(0)
    if opt == "--max-subwalk-size":
        max_subwalk_size = int(sys.argv.pop(0))
        assert max_subwalk_size >= 64
        #print(f"max_subwalk_size: {max_subwalk_size}")
    else:
        # can't reach
        print(f"error: unsupported option {opt}")
        raise AssertionError()
os.environ["ASM"] = "/Users/subarutaro/work/pfcomp/gpfn2/assembler/asm"
os.environ["PACKER"] = "/Users/subarutaro/work/pfcomp/gpfn2/emulator/build/grape_pfn/packer/packer_main"

fl72 = lambda v:hex(struct.unpack('>Q', struct.pack('>d', v))[0]).split("x")[1]+"00"

with open("codehead.json","r") as f:
    jsonstr = f.read()
    load_obj = json.loads(jsonstr)
    
    regs_dic = load_obj["regs_dic"]
    core_instlist_head = load_obj["core_instlist_head"]
    imm_instlist_head = load_obj["imm_instlist_head"]
    decllist = load_obj["decllist"]
#print(regs_dic)
ChI = len( regs_dic['I'] )
ChJ = len( regs_dic['J'] )
ChF = len( regs_dic['F'] )

max_element_size=16
min_element_size=64
io_size= [0,0,0]
var_type_list=[[],[],[]]
for index, io in enumerate(['I','J','F']):
    for v in regs_dic[io]:
        var, tp = v.split(':')
        size=int(tp[1:])
        io_size[index] += size
        var_type_list[index].append([var,tp])
        min_element_size=min(min_element_size,size)
        max_element_size=max(max_element_size,size)

max_inst_prec = 16
min_inst_prec = 64
for insts in core_instlist_head:
    for op, _, _, _ in insts:
        op_code = op.split(':')[0]
        if op_code in ['nop','noforward','mask','mend','']:
            continue
        tp = op.split(':')[1]
        size = int(tp[1:])
        max_inst_prec = max(max_inst_prec,size)
        min_inst_prec = min(min_inst_prec,size)

is_multi_prec = (max_inst_prec != min_inst_prec)


#NI, NJ= 256, 4
#NI, NJ = 512, 2 # NIpL1B, NJpPE
simd_width=64//min_inst_prec
ni_per_loop = simd_width*cycle_per_step
nj_per_loop = simd_width
NI = ((LMEM_SIZE-LMEM_BANK_OFFSET) // (max(io_size[0],io_size[2])//SW_BIT_SIZE) // 64) * num_pe_per_l1b
if max_subwalk_size > 0:
    NI = min( NI, simd_width*max_subwalk_size )
NI = (NI//(simd_width*l2bm_chunk_size))*(simd_width*l2bm_chunk_size) # NI must be multiple of simd_width*l2bm_chunk_size
#NJ = 2 if min_element_size == 64 else 4
NJ = ((min((GRF_SIZE-GRF_BANK_OFFSET) // (SW_IN_LW*(io_size[1]//SW_BIT_SIZE)), DRAM_OFFSET_SIZE // (SW_IN_LW*(io_size[1]//SW_BIT_SIZE)*num_pe)))//simd_width)*simd_width

#print("ChX:",ChI,ChJ,ChF)
#print("(NI,NJ):",NI,NJ)

EPI_OFFSET   = [ 0 for i in range(ChI+1) ]
EPJ_OFFSET   = [ 0 for i in range(ChJ+1) ]
FORCE_OFFSET = [ 0 for i in range(ChF+1) ]
EPI_OFFSET[0]   = LMEM_BANK_OFFSET
EPJ_OFFSET[0]   = GRF_BANK_OFFSET
FORCE_OFFSET[0] = LMEM_BANK_OFFSET
DRAM_EPI_OFFSET   = [ 0 for i in range(ChI+1) ]
DRAM_EPJ_OFFSET   = [ 0 for i in range(ChJ+1) ]
DRAM_FORCE_OFFSET = [ 0 for i in range(ChF+1) ]
DRAM_EPI_OFFSET[0] = 0
DRAM_EPJ_OFFSET[0] = DRAM_IN_J_OFFSET
DRAM_FORCE_OFFSET[0] = DRAM_OUT_OFFSET
for i in range(len(var_type_list[0])):
    var=var_type_list[0][i][0]
    tp=var_type_list[0][i][1]
    EPI_OFFSET[i+1] = EPI_OFFSET[i] + NI*(int(tp[1:3])//32)
    DRAM_EPI_OFFSET[i+1] = DRAM_EPI_OFFSET[i] + NI//(64//int(tp[1:3]))*num_l1b
    if(EPI_OFFSET[i+1]%2 != 0):
        EPI_OFFSET[i+1] += 1
for i in range(len(var_type_list[1])):
    var=var_type_list[1][i][0]
    tp=var_type_list[1][i][1]
    EPJ_OFFSET[i+1] = EPJ_OFFSET[i] + NJ*(int(tp[1:3])//32)
    DRAM_EPJ_OFFSET[i+1] = DRAM_EPJ_OFFSET[i] + num_l1b * (NJ*num_pe_per_mab*num_mab_per_l1b//(64//int(tp[1:3])))
    if(EPJ_OFFSET[i+1]%2 != 0):
        EPJ_OFFSET[i+1] += 1
for i in range(len(var_type_list[2])):
    var=var_type_list[2][i][0]
    tp=var_type_list[2][i][1]
    FORCE_OFFSET[i+1] = FORCE_OFFSET[i] + NI*(int(tp[1:3])//32)
    DRAM_FORCE_OFFSET[i+1] = DRAM_FORCE_OFFSET[i] + NI//(64//int(tp[1:3]))*num_l1b
    if(FORCE_OFFSET[i+1]%2 != 0):
        FORCE_OFFSET[i+1] += 1
#print(EPI_OFFSET)
#print(EPJ_OFFSET)
#print(FORCE_OFFSET)
#print(DRAM_EPJ_OFFSET)

assert EPI_OFFSET[ChI]   <= LMEM_SIZE
assert EPJ_OFFSET[ChJ]   <= GRF_SIZE
assert FORCE_OFFSET[ChF] <= LMEM_SIZE

raw_op_converter = {
    "add:F64":"dvadd",
    "sub:F64":"dvadd",
    "fmau:F64":"dvfmau",
    "fmad:F64":"dvfmad",
    "set:F64":"lpassa",
    "rsqrt:F64":"drsqrt",
    "max:F64":"dmax",
    "min:F64":"dmin",
    "mulu:F64":"dvmulu",
    "muld:F64":"dvmuld",
    "fmsu:F64":"dvfmau",
    "fmsd:F64":"dvfmad",
    "fmiu:F64":"dvfmau",
    "fmid:F64":"dvfmad",
    "add:F32":"fvadd",
    "sub:F32":"fvadd",
    "fma:F32":"fvfma",
    "set:F32":"lpassa",
    "rsqrt:F32":"frsqrt",
    "max:F32":"fmax",
    "min:F32":"fmin",
    "mul:F32":"fvmul",
    "fms:F32":"fvfma",
    "fmi:F32":"fvfma",
    "add:S32":"add",
    "sub:S32":"sub",
    "add:U32":"add",
    "sub:U32":"sub",
    "nop:":"nop",
    "noforward:":"noforward",
    "lnot:BOOL":"ilnot",
    "land:BOOL":"land",
    "and:F64":"land",
    "and:S64":"land",
    "and:F32":"iand",
    "and:S32":"iand",
    "xor:F64":"lxor",
    "xor:S64":"lxor",
    "xor:F32":"ixor",
    "xor:S32":"ixor",
    "dec:F64":"ldec",
    "dec:S64":"ldec",
    "dec:F32":"idec",
    "dec:S32":"idec",
    "rotate:F32":"lbsr",
    "rotate:F16":"lbsr",
    "imm:F64":"imm",
    "imm:F32":"imm",
    "imm:F16":"imm",
}
mau_inst = [
    "add",
    "sub",
    "mul",
    "mulu",
    "muld",
    "fma",
    "fmau",
    "fmad",
    "fms",
    "fmsu",
    "fmsd",
    "fmi",
    "fmiu",
    "fmid",
        ]

# generate dict for memory address
id_to_reg = {}
for regtype in "IJFRSMNO":
    addr = 0
    for insts in core_instlist_head:
        for inst in insts:
            op, reads, writes, _ = inst
            for reg in reads+writes:
                if reg in id_to_reg.keys():
                    continue
                if regtype == reg.split('.')[0]:
                    items = reg.split(':')
                    regid = items[0].split('.')[1]
                    tp = items[1]
                    simd = int(items[2]) if len(items)>2 else None
                    size = int(tp[1:]) if tp != "BOOL" else 0
                    s_or_l = "long"
                    if tp == "F32" and simd != None:
                        s_or_l = "short" 
                    elif op.split(':')[1] == "F32" and tp == "F64":
                        s_or_l = "double_long"
                        assert False
                    if regtype in "IJFO":
                        id_to_reg[reg] = [s_or_l,int(regid),simd]
                    else:
                        reg_wo_simd = regtype + "." + regid + ":" + tp
                        if reg_wo_simd not in id_to_reg.keys():
                            id_to_reg[reg_wo_simd] = ["long" if s_or_l == "short" else s_or_l ,addr,None] # non-simd var can't be short
                            addr += 2*cycle_per_step
                        if simd != None:
                            _, base, _ = id_to_reg[reg_wo_simd]
                            if tp == "F32":
                                id_to_reg[reg] = [s_or_l,base+simd,simd]
                            else:
                                id_to_reg[reg] = [s_or_l,base+2*simd,simd]

def getEpiRegAddress(i,c):
    tp = var_type_list[0][c][1]
    word_length = 2 if '64' in tp else 1
    addr = EPI_OFFSET[c] + word_length*i
    return addr

def getEpjRegAddress(j,c,ISX,tp='F64',s=None):
    tp = var_type_list[1][c][1]
    simd = 2 if '32' in tp else 1
    assert NJ%simd == 0
    BANKOFFSET = (0 if ISX else EPJ_OFFSET[ChJ] - GRF_BANK_OFFSET)
    #addr = BANKOFFSET+2*((ChJ*(j//simd))+c)
    addr = BANKOFFSET+EPJ_OFFSET[c] + 2*(j//simd)
    if s != None and simd == 1:
        if j%2 == 0:
            #addr = BANKOFFSET+2*((ChJ*((j+s)//simd))+c)
            addr = BANKOFFSET+EPJ_OFFSET[c] + 2*((j+s)//simd)
        else:
            #addr = BANKOFFSET+2*((ChJ*((j-s)//simd))+c)
            addr = BANKOFFSET+EPJ_OFFSET[c] + 2*((j-s)//simd)

    #print(addr,BANKOFFSET,GRF_BANK_OFFSET,NJ,ChJ,j,c,ISX)
    assert addr <= GRF_SIZE - 2
    assert addr%2 == 0
    return addr

def getForceRegAddress(i,c):
    tp = var_type_list[2][c][1]
    word_length = 2 if '64' in tp else 1
    addr = FORCE_OFFSET[c] + word_length*i
    return addr

def getForceTmpRegAddress(regid,i=None):
    if i != None:
        assert i < simd_width
    addr = GRF_BANK_OFFSET + 8*simd_width*regid
    if i != None:
        addr += 8*i
    assert (addr + 2*cycle_per_step)<= GRF_SIZE
    return addr

def getForceTmpRegName(regid,i=0):
    reg = None
    addr = getForceTmpRegAddress(regid,i)
    reg = f"$ls{addr}v2"
    return reg

def getForceRegName(regid,idx):
    addr = getForceRegAddress(idx,regid)
    assert (addr + 2*cycle_per_step) <= LMEM_SIZE
    reg = f"$ln{addr}v2"
    return reg 

def getRegName(regtype, addr, i_index, j_index, ISX, NI, ChJ, ChF, tp_, s_or_l = "long", simd = None):
    reg = None
    GRF_SIZE_for_f = ChF*2*cycle_per_step
    prefix = "$"
    if s_or_l == "long":
        prefix += "l"
    elif s_or_l == "double_long":
        prefix += "ll"
    suffix = "v2"
    if s_or_l == "short":
        if simd != None:
            suffix = "v2"
        else:
            suffix = "v1"
    elif (s_or_l == "long" and simd != None) or s_or_l == "double_long":
        suffix = "v4"
    if regtype == "R":
        #print(addr,GRF_BANK_OFFSET,GRF_SIZE_for_f)
        assert addr <= GRF_BANK_OFFSET - 2*cycle_per_step
        assert addr <= GRF_SIZE - 2*cycle_per_step
        reg = f"{prefix}r{addr}{suffix}"
    if regtype == "S":
        assert addr <= GRF_BANK_OFFSET - 2*cycle_per_step
        assert addr <= GRF_SIZE - 2*cycle_per_step
        reg = f"{prefix}s{addr}{suffix}"
    if regtype == "M":
        assert addr <= LMEM_BANK_OFFSET - 2*cycle_per_step
        assert addr <= LMEM_SIZE - 2*cycle_per_step
        reg = f"{prefix}m{addr}{suffix}"
    if regtype == "N":
        assert addr <= LMEM_BANK_OFFSET - 2*cycle_per_step
        assert addr <= LMEM_SIZE - 2*cycle_per_step
        reg = f"{prefix}n{addr}{suffix}"
    if regtype == "I":
        #BANKOFFSET = LMEM_BANK_OFFSET
        #reg = f"$lm{2*128+2*(i_index*ChI)+2*regid}v{2*ChI}" # original
        #addr = BANKOFFSET + 2*(NI*regid + i_index)
        if simd != None:
            i_index += simd
        addr = getEpiRegAddress(i_index,addr)
        assert addr <= LMEM_SIZE - 2*cycle_per_step
        reg = f"{prefix}m{addr}{suffix}"
    if regtype == "J":
        addr = getEpjRegAddress(j_index,addr,ISX,tp_,simd)
        assert addr <= GRF_SIZE - 2
        reg = f"{prefix}r{addr}v0"
    if regtype == "F":
        addr = getForceTmpRegAddress(addr,simd)
        if simd != None:
            addr -= 6*simd
        reg = f"{prefix}s{addr}{suffix}"
    if regtype == "A":
        reg = "$aluf"
    if regtype == "U":
        reg = "$mauf"
    if regtype == "Z":
        reg = "$nowrite"
    if regtype == "O":
        reg = f"$omr{addr}"
    if reg == None:
        print(regtype, addr)
        assert False
    return reg

isImm = lambda x: x.split('#')[0] in ["VF","VI"]
isConstant = lambda x: x.split('#')[0] in ["C"]
def immediate_value(reg):
    _, tmp = reg.split('#')
    val, tp = tmp.split(':')
    if val[-1] == 'f':
        val = val.split('f')[0]
    assert tp in ["F32","S32","U32"]
    prefix = 'f'
    if tp == "U32":
        assert int(val) >= 0 and int(val) <= 4294967295
        prefix = 'u'
    elif tp == "S32":
        assert int(val) >= -2147483648 and int(val) <= 2147483647
        prefix = "i"
    return f"{prefix}\"{val}\""

def constant_value(reg):
    val = reg.split('#')[1].split(':')[0]
    if val == "MSB":
        return "$msb1"
    #elif val == "PEID":
    #    return "$peid"
    #elif val == "SUBPEID":
    #    return "$subpeid"
    else:
        # can't reach
        print(reg)
        raise AssertionError()

forwardings=["A","U","Z"]
def gen_Main(core_instlist_head, regs_dic, NI, NJ, ISX=True): # NI:NI per L1B, NJ:NJ per PE
    #assert NI==256 and NJ == 4

    codeM = Code()
    for i_index in range(0,NI,ni_per_loop):
        for regid in range(ChF):
            # this passa can be merged to other MAU instruction?
            # pass F on LN to LS
            tp = regs_dic['F'][regid].split(':')[1]
            if simd_width > 1 and int(tp[1:]) == 64:
                for i in range(simd_width):
                    src=getForceRegName(regid,i_index+i*cycle_per_step)
                    dst=getForceTmpRegName(regid,i)
                    codeM[  f"dpassa {src} {dst}"  ]
            else:
                src=getForceRegName(regid,i_index)
                dst=getForceTmpRegName(regid)
                codeM[  f"dpassa {src} {dst}"  ]
        
        for j_index in range(NJ):
            
            #CORE PART START
            for insts in core_instlist_head:
                vsm = ""
                for op, reads, writes, md in insts:
                    if op == "":
                        continue
                    op_ = op.split(":")[0]
                    tp_ = op.split(":")[1]
                    regs = reads + writes
                    raw_regs = []
                    for rid, reg in enumerate(regs):
                        if isImm(reg):
                            reg = immediate_value(reg)
                        elif isConstant(reg):
                            reg = constant_value(reg)
                        else:
                            regtype = get_reg_name(reg)
                            tp = get_reg_type(reg)
                            s_or_l, addr, simd = id_to_reg.get(reg,["long",0,0])
                            reg = getRegName(regtype, addr, i_index, j_index, ISX, NI, ChJ, ChF, tp_, s_or_l, simd)
                            mask = f"/$imr{md}" if md>0 and reg != "$nowrite" else ""
                            if regtype not in ["A","U","Z"] and tp_ == "F64" and tp == "F32" and rid < len(reads) and op_ in mau_inst:
                                reg += 'e'
                            if s_or_l == "double_long":
                                reg += 'r'
                        vlen = 4
                        reg = reg+mask if rid >= len(reads) else reg
                        if op_ in ["fmiu","fmid","fmi"]: # -a*b + c
                            reg = "-"+reg if rid == 0 else reg
                        if op_ in ["sub"]:
                            reg = "-"+reg if rid == 1 else reg
                        if op_ in ["fmsu","fmsd","fms"]: # a*b - c
                            reg = "-"+reg if rid == 2 else reg
                        raw_regs.append(reg)

                    #print ( [op, reads, writes, md] )
                    sole = False if md==0 else True

                    if op in raw_op_converter:
                        raw_op = raw_op_converter[op]
                        if vsm != "":
                            vsm += "; "
                        opt = op.split(":")[1]
                        if len(writes)>0:
                            wtp0 = writes[0].split(":")[1]
                            assert all([ w.split(":")[1] == wtp0 or w.split(".")[0] == "O" for w in writes])
                            if wtp0 == "F32" and opt == "F64":
                                if raw_op in ["dvadd","dvmulu","dvmuld", "dvfmau", "dvfmad"]:
                                    raw_op += "r"

                        vsm += " ".join([raw_op,*raw_regs])
                        #print (" →", codeM.ops[-1] )
                    #elif op == "cmp_mw:F64":
                    #    raw_op = "dvadda"
                    #    codeM[   " ".join(  [raw_op,*raw_regs] + [f"$omr{1+md}"]  ) ,mask,vlen,sole  ]
                    #    #print (" →", codeM.ops[-1] )
                    #elif op == "cmp_fl:BOOL":
                    #    raw_op = "imm"
                    #    codeM[   " ".join(  [raw_op, 'i"0"' ,*raw_regs]  ) ,mask,vlen,sole  ]
                    #    #print (" →", codeM.ops[-1] )
                    #elif op == "cmp_tr:BOOL":
                    #    raw_op = "imm"
                    #    codeM[   " ".join(  [raw_op, 'i"1"' ,*raw_regs]  ) ,mask,vlen,sole  ]
                    #    #print (" →", codeM.ops[-1] )
                    #elif op == "mwrt:BOOL":
                    #    raw_op = "idec"
                    #    codeM[   " ".join(  [raw_op,*raw_regs] + [f"$omr{1+md}"]  ) ,mask,vlen,sole  ]
                    #     #print (" →", codeM.ops[-1] )
                    else:
                        print (" →", "\033[91m {}\033[00m" .format("Not Implemented Error"), op )
                        assert 1==0
                codeM[vsm]
            #CORE PART END
        
        #codeM[ "nop" ]
        for regid in range(ChF):
            tp = regs_dic['F'][regid].split(':')[1]
            if simd_width > 1 and int(tp[1:]) == 64:
                for i in range(simd_width):
                    src=getForceTmpRegName(regid,i)
                    dst=getForceRegName(regid,i_index+i*cycle_per_step)
                    codeM[  f"dpassa {src} {dst}"  ]
            else:
                src=getForceTmpRegName(regid)
                dst=getForceRegName(regid,i_index)
                codeM[  f"dpassa {src} {dst}"  ]
        codeM[ "nop" ]
        codeM[ "nop" ]
    return codeM

def DRAM2LOCAL_EPJ(code,dram_offset,dram_increment_per_l2b,localmem_offset,target_localmem="n"):
    #Send an element per PE from DRAM to localmem. EPJ data is not shared at all.
    #This operation send 1LW per PE, 4LW per MAB, 4*16LW per L1B, 4*16*8LW per L2B, 4*16*8*8LW per Device.
    #This requires to set 'target_localmem' as 'm' or 'n'
    for group_id in range(num_group):
        for l2b_id in range(num_l2b_per_group):
            index = num_l2b_per_group * group_id + l2b_id
            tag = "i" + format(index+1,'02x')
            code[f'mvp/n{dram_increment_per_l2b}{tag} $d{dram_offset + dram_increment_per_l2b*(index)}@0 $lc0@{group_id}.{l2b_id}' ]
            code[f'nop; wait {tag}']

    [ code[f"l2bmb@{l1b_id} $lc{l2bm_chunk_size*l1b_id} $lb{0}"] for l1b_id in range(num_l1b_per_l2b) ]
    [ code['nop'] for i in range(4)]
    code[f"l1bmd $lb{0} $l{target_localmem}{localmem_offset}/1000"]

def gen_Sub( regs_dic, NI, NJ, ISX=True, dstart=0):
    offset = dstart + (0 if not ISX else DRAM_OFFSET_SIZE)

    codeS = Code()
    for c in range(ChJ):
        tp = var_type_list[1][c][1]
        size = int(tp[1:])
        interval = 64//size
        for j in range(0,NJ,interval):
            lmem_offset = getEpjRegAddress(j,c,not ISX)
            doffset = offset + DRAM_EPJ_OFFSET[c] + num_pe*j//interval
            # data on dram are in order below:
            # nj[0][0], ... , nj[0][NJ*num_pe], nj[1][0], ... , nj[1][NJ*num_pe], ... , nj[ChJ-1][0], ...
            DRAM2LOCAL_EPJ(codeS,dram_offset=doffset, dram_increment_per_l2b= num_pe_per_l2b, localmem_offset= lmem_offset, target_localmem="r")
    return codeS

def overlap_pe_alu_inst(code):
    assert False # under construction
    codeA = Code()

    alu_insts = ["rsqrt","max","min","floor","bsr","bsl"]
    op_list = list(reversed(copy.deepcopy(code.ops))) # !REVERSED

    for i in range(100000):
        if len(op_list) == 0:
            break
        op = op_list.pop()
        if ';' in op[0]:
            codeA.ops.append(op) # pe + alu inst
            continue
        if any( [(alu_inst in op[0]) for alu_inst in alu_insts] ): # alu inst only
            alu_ops = op[0].split()
            alu_write = alu_ops.pop()
            alu_reads = []
            while len(alu_ops) > 1:
                alu_reads.append(alu_ops.pop())
            isOverlapped = False
            opsTmp = []
            while (not isOverlapped) and (len(codeA.ops)>0):
                ops = codeA.ops.pop()
                write = ops[0].split()[-1]
                if write in alu_reads:
                    codeA.ops.append(ops)
                    break
                opsTmp.append(ops)
            while (not isOverlapped) and (len(opsTmp)>0):
                op_bufferd = opsTmp.pop()
                ops = op_bufferd.split()
                if not isOverlapped:
                    write = ops.pop()
                    reads = []
                    while len(alu_ops) > 1:
                        reads.append(ops.pop())
                    
                    
                codeA.ops.append(op_buffered)
            if not isOverlapped:
                codeA.ops.append(op)
            continue
        # pe inst only
        codeA.ops.append(op)
    return codeA

def Merge_MS ( codeM, codeS ):
    codeA = Code()
    
    def IsMVInst (op):
        o,m,v,s = op
        return o.strip()[:2] == 'mv'
    
    opM_list = list(reversed(copy.deepcopy(codeM.ops))) # !REVERSED
    opS_list = list(reversed(copy.deepcopy(codeS.ops))) # !REVERSED

    single_operand_insts = [
            "imm",
            "msl", "msr",
            "ipassa", "lpassa",
            "iinc", "linc",
            "idec", "ldec",
            "inot", "lnot",
            "ilnot", "llnot",
            "frsqrt", "drsqrt", "hrsqrt",
            "ffloor", "dfloor", "hfloor",
            "fftoi", "dftoi", "hftoi",
            ]
    triple_operand_insts = [
            "fvfma", "dvfmau", "dvfmad"
            ]

    mv_timer = 0
    for i in range(1000000):
        mv_timer -= 1
        if i==999999: 
            assert 1==0 # something wrong leads to infinity loop!
        if len(opM_list)==0 or len(opS_list)==0: #either list is empty!
            break
        opM = opM_list.pop()
        opS = opS_list.pop()
        doWait = True if ('wait' in opS[0]) and (mv_timer <= 0) else False
        waitOp = opS[0].split(';')[1] if doWait else ''
        if doWait:
            opM[0] += ';' + waitOp
        if IsMVInst(opS):
            codeA.ops.append(opS)
            mv_timer = 64 + (int(opS[0].split('n')[1].split('i')[0])//64) # assuming 16 LW/cycle * 4 cycle/step
            opM_list.append(opM)
            continue
        if opM[0] == 'nop':
            if ('wait' in opS[0]) and (mv_timer > 0):
                codeA.ops.append(opM)
                opS_list.append(opS)
            else:
                codeA.ops.append(opS)
                assert not doWait
            continue
        else:
            if opS[0] == 'nop':
                codeA.ops.append(opM)
                assert not doWait
                continue
            if ('l2bm' in opS[0]):
                o,m,v,s = opM
                codeA.ops.append([opS[0]+";"+o,m,v,s])
                continue
            if ('l1bm' in opS[0]):
                o,m,v,s = opM
                insts = o.split(';')
                canOverlap = True
                for inst in insts:
                    op_ops = inst.split()
                    op = op_ops[0]
                    ops = op_ops[1:]
                    num_reads = 2
                    if op in single_operand_insts:
                        num_reads = 1
                    elif op in triple_operand_insts:
                        num_reads = 3
                    reads = ops[:num_reads]
                    writes = ops[num_reads:]
                    canOverlap = canOverlap and all( [ not ('$lr' in reg or '$r' in reg) for reg in writes])
                    canOverlap = canOverlap and all( ['/$imr' not in reg for reg in writes])
                if canOverlap:
                    codeA.ops.append([opS[0]+"; "+o,m,v,s])
                else:
                    codeA.ops.append(opM)
                    opS_list.append(opS)
                continue
            else:
                codeA.ops.append(opM)
                if not doWait:
                    opS_list.append(opS)
                continue
            assert 1==0
    #print ( len(opM_list), len(opS_list), ":Merge remain (Main, Sub)" )
    for i in range(len(opM_list)):
        codeA.ops.append(opM_list.pop())
    for i in range(len(opS_list)):
        codeA.ops.append(opS_list.pop())
    return codeA 

def cat_MS(codeM, codeS):
    codeA = Code()
    opM_list = list(reversed(copy.deepcopy(codeM.ops))) # !REVERSED
    opS_list = list(reversed(copy.deepcopy(codeS.ops))) # !REVERSED
    for i in range(len(codeM.ops)):
        codeA.ops.append(opM_list.pop())
    for i in range(len(opS_list)):
        codeA.ops.append(opS_list.pop())
    return codeA


def DRAM2LOCAL_EPI(code,dram_offset,dram_increment_per_l2b,localmem_offset, target_localmem):
    #Send from DRAM to localmem. EPI data is shared in L1B.
    #This operation send 64LW per PE, 64LW per L1B, 64*8LW per L2B, 64*8*16LW per Device.
    #also, both localmen LM,LN is shared.
    for group_id in range(num_group):
        for l2b_id in range(num_l2b_per_group):    
            index = num_l2b_per_group*group_id + l2b_id
            tag = "i" + format(index+1,'02x')
            code[f'mvp/n{dram_increment_per_l2b}{tag} $d{dram_offset + dram_increment_per_l2b*index}@0 $lc0@{group_id}.{l2b_id}' ]
            code[f'nop; wait {tag}']
    [ code[f"l2bmb@{l1b_id} $lc{l2bm_chunk_size*l1b_id} $lb{0}"] for l1b_id in range(num_l1b_per_l2b) ]
    [ code['nop'] for i in range(4)]
    [ code[f"l1bmp $lb{4*i} $l{target_localmem}{localmem_offset+2*(4*i)}v2"] for i in range(l2bm_chunk_size//cycle_per_step)]

def LOCAL2DRAM_FORCE(code,dram_offset,dram_increment_per_l2b, target_localmem, element_size, c, simd, nloop):
    assert target_localmem in ['m','n']
    interval = 2*cycle_per_step
    for j in range(nloop) :
        localmem_offset = getForceRegAddress(simd*l2bm_chunk_size*j, c)
        for i in range( l2bm_chunk_size//cycle_per_step ):
            Foff = localmem_offset + interval*i
            F = f'l{target_localmem}{Foff}'
            # accumulate among PEs
            if 'F64' in element_size:
                code[f'                             msl ${F}v2 $nowrite','',4,True]
                code[f'dvadd $aluf ${F}v2 $nowrite; msl $aluf  $nowrite','',4,True]
                code[f'dvadd $aluf $mauf  $nowrite; msl $aluf  $nowrite','',4,True]
                code[f'dvadd $aluf $mauf  $nowrite                     ','',4,True]
                # write accumulated to each PE's Foff adress
                code[f'dmwrite $mauf $lx0','',4,True]
                code[f'dmread $lx0 ${F}'  ,'',4,True]
            elif 'F32' in element_size:
                code[f'                             msl ${F}v2 $nowrite','',4,True]
                code[f'fvadd $aluf ${F}v2 $nowrite; msl $aluf  $nowrite','',4,True]
                code[f'fvadd $aluf $mauf  $nowrite; msl $aluf  $nowrite','',4,True]
                code[f'fvadd $aluf $mauf  $nowrite                     ','',4,True]
                # write accumulated to each PE's Foff adress
                code[f'dmwrite $mauf $lx0','',4,True]
                code[f'dmread $lx0 ${F}'  ,'',4,True]
            else:
                print(f"LOCAL2DRAM_FORCE does not support #{element_size}")
                assert False
        
            [ code['nop'] for i in range(2) ]    #################

        for i in range(4):
            addr = localmem_offset+interval*cycle_per_step*i
            if 'F64' in element_size:
                code[f'l1bmrdfadd $l{target_localmem}{addr}v{interval} $lb{16*i + l2bm_chunk_size*j}']
            elif 'F32' in element_size:
                code[f'l1bmrffadd $l{target_localmem}{addr}v{interval} $lb{16*i + l2bm_chunk_size*j}']
            else:
                print(f"LOCAL2DRAM_FORCE does not support #{element_size}")
                assert False
    if c > 0:
        for group_id in range(num_group):
            for l2b_id in range(num_l2b_per_group):
                index = num_l2b_per_group*group_id + l2b_id
                tag = "i" + format(index+1,'02x')
                code[f'nop; wait {tag}']
    else:
        [ code['nop'] for i in range(4)]

    for j in range(nloop) :
        [ code[f"l2bm@{l1b_id} $lb{l2bm_chunk_size*j} $lc{nloop*l2bm_chunk_size*l1b_id + l2bm_chunk_size*j}"] for l1b_id in range(8) ] #L2BM write 64x4x8

    for group_id in range(num_group):
        for l2b_id in range(num_l2b_per_group):
            index = num_l2b_per_group*group_id + l2b_id
            tag = "i" + format(index+1,'02x')
            code[f'mvp/n{nloop*l2bm_chunk_size*num_l1b_per_l2b}{tag} $lc0@{group_id}.{l2b_id} $d{dram_offset+dram_increment_per_l2b*index}@0 ' ]
            if c == ChF-1:
                code[f'nop; wait {tag}']

def double_to_hex(f):
    return hex(struct.unpack('<Q', struct.pack('<d', f))[0])

def gen_Head(regs_dic, NI, NJ, offset):
    
    codeA = Code()
    
    for c in range(ChI) :
        codeA[f"# copy {var_type_list[0][c][0]} (EPI{c}) "]
        tp = var_type_list[0][c][1]
        size = int(tp[1:])
        nloop = NI // ((64//size)*l2bm_chunk_size)
        for j in range(nloop):
            # data of ni is in order of :
            # ni[0][0], ... , ni[0][NI*num_l1b-1], ni[1][0], ..., ni[ChI-1][0], ... , ni[ChI-1][NI*num_l1b-1]
            #doffset = DRAM_IN_OFFSET+((NI//(64//size))*num_l1b*c + num_l1b*64*j)
            doffset = offset + DRAM_EPI_OFFSET[c] + num_l1b*64*j
            lmem_offset = getEpiRegAddress(((64//size)*64)*j,c)
            #assert doffset < DRAM_OUT_OFFSET
            codeA[f"#DRAM2LOCAL_EPI c={c} j={j} dram_offset={doffset} lmem_offset={lmem_offset}"]
            DRAM2LOCAL_EPI(codeA, dram_offset=doffset, dram_increment_per_l2b=num_l1b_per_l2b*l2bm_chunk_size, localmem_offset=lmem_offset, target_localmem="m")
    [codeA["nop"] for i in range(4)]
    #print ( "@", len(codeA.ops) )
    
    #NF clear
    for c in range( ChF ):
        codeA[f"# clear {var_type_list[2][c][0]} (FORCE{c}) "]
        tp = var_type_list[2][c][1]
        size = int(tp[1:])
        interval = (64//size)*cycle_per_step
        for i in range( 0, NI, interval ):
            dst=getForceRegName(c,i)
            codeA[f'zero {dst}']
    #print ( "@", len(codeA.ops) )
    #NJ clear ( All Clear of Bank 0 and Bank 1 )
    for k in range( (GRF_SIZE - GRF_BANK_OFFSET) // (2*cycle_per_step) ):
        codeA[f'zero $lr{GRF_BANK_OFFSET + 2*cycle_per_step*k}v2']
    #print ( "@", len(codeA.ops) )
        
    #imm (not repeat required)
    [codeA["nop"] for i in range(4)]
    for reg, val in imm_instlist_head:
        regtype = get_reg_name(reg)
        regid = int(get_reg_id(reg)) if regtype != "T" else None
        valtype = val.split(':')[1]
        loc = regtype.lower()
        _, addr, _ = id_to_reg[reg]
        if valtype == 'F64':
            v = float(val.split(":")[0])
            h = int(double_to_hex(v)[0:10],0)
            #print(v,h)
            codeA[f'immu i"{h}" $nowrite']
            codeA[f'lpassa $aluf $l{loc}{addr}v2']
        elif valtype == 'F32':
            v = val.split(":")[0].split("f")[0]
            codeA[f'imm f"{v}" $nowrite']
            codeA[f'lpassa $aluf $l{loc}{addr}v2']
        elif valtype == 'U64':
            v = int(val.split(":")[0])
            assert v < 2147483648 # 2**31
            codeA[f'imm i"{v}" $nowrite']
            codeA[f'lpassa $aluf $l{loc}{addr}v2']
            codeA[f'zero   ${loc}{addr}v2']
        elif valtype == 'U32':
            v = int(val.split(":")[0])
            assert v < 2147483648 # 2**31
            codeA[f'imm i"{v}" $nowrite']
            codeA[f'lpassa $aluf $l{loc}{addr}v2']
        else:
            assert False
            
    [codeA["nop"] for i in range(4)]
    #print ( "@", len(codeA.ops) )

    # first NJ transfer
    for c in range(ChJ):
        codeA[f"# copy {var_type_list[1][c][0]} (EPJ{c}) "]
        tp = var_type_list[1][c][1]
        size = int(tp[1:])
        interval = 64//size
        for j in range(0,NJ,interval):
            lmem_offset = getEpjRegAddress(j,c,True)
            # data on dram are in order below:
            # nj[0][0], ... , nj[0][NJ*num_pe], nj[1][0], ... , nj[1][NJ*num_pe], ... , nj[ChJ-1][0], ...
            doffset = offset + DRAM_EPJ_OFFSET[c] + num_pe*j//interval
            assert doffset < offset + DRAM_IN_J_OFFSET + DRAM_OFFSET_SIZE
            codeA[f"#DRAM2LOCAL_EPJ c={c} j={j} dram_offset={doffset} lmem_offset={lmem_offset}"]
            DRAM2LOCAL_EPJ(codeA,dram_offset=doffset, dram_increment_per_l2b= num_pe_per_l2b, localmem_offset= lmem_offset,target_localmem="r")

    return codeA

def gen_Tail(regs_dic, NI, NJ, offset):
    codeA = Code()
    
    for c in range( ChF ) :
        codeA[f"# copy {var_type_list[2][c][0]} (FORCE{c})"]
        tp = var_type_list[2][c][1]
        size = int(tp[1:])
        simd = 64 // size
        nloop = NI // ((64//size)*l2bm_chunk_size)
        doffset = offset + DRAM_FORCE_OFFSET[c]
        codeA[f"#LOCAL2DRAM_FORCE c={c} dram_offset={doffset}"]
        LOCAL2DRAM_FORCE(codeA, dram_offset=doffset, dram_increment_per_l2b=num_l1b_per_l2b*(NI//simd),target_localmem="n", element_size=var_type_list[2][c][1], c=c, simd=simd, nloop=nloop)
    #print ( "@", len(codeA.ops) )

    return codeA

for i in range(min(DRAM_SIZE//DRAM_KERNEL_OFFSET,MAX_KERNEL_LAUNCH)):
    offset = i*DRAM_KERNEL_OFFSET
    ISX = True
    codeM = gen_Main(core_instlist_head, regs_dic, NI=NI, NJ= NJ, ISX=ISX)
    codeS = gen_Sub( regs_dic, NI=NI, NJ= NJ, ISX=ISX, dstart=offset)
    codeMS = Merge_MS(codeM,codeS)
    #codeMS = cat_MS(codeM,codeS)

    ISX = False
    codeMr = gen_Main(core_instlist_head, regs_dic, NI=NI, NJ= NJ, ISX=ISX)
    codeSr = gen_Sub( regs_dic, NI=NI, NJ= NJ, ISX=ISX, dstart=offset)
    codeMSr = Merge_MS(codeMr,codeSr)
    #codeMSr = cat_MS(codeMr,codeSr)

    name = f"packs/nbY1_MS{i}"
    WriteTextCode(  ToTextCode(codeMS) , name=name )

    name = f"packs/nbY1_MSR{i}"
    WriteTextCode(  ToTextCode(codeMSr) , name=name )

    #print ( len(codeM.ops), len(codeS.ops), len(codeMS.ops) )

    codeH = gen_Head( regs_dic, NI, NJ, offset)
    name = f"packs/nbY1_H{i}"
    WriteTextCode(  ToTextCode(codeH)  , name=name )
    #print ( len(codeH.ops) )

    codeT = gen_Tail( regs_dic, NI, NJ, offset)
    name = f"packs/nbY1_MT{i}"
    WriteTextCode(  ToTextCode(codeM)+ToTextCode(codeT)  , name=name )
    name = f"packs/nbY1_MRT{i}"
    WriteTextCode(  ToTextCode(codeMr)+ToTextCode(codeT)  , name=name )
    #print ( len(codeT.ops) )
print("finished")
