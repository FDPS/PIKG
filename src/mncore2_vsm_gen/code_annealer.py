from simanneal import Annealer
import random

from generate_core_instlist import get_reg_id, get_reg_name, get_reg_type

def code_annealer(core_instlist_gen, Rnum, Ilen, nloops=5):
    res_li = []
    for i in range(nloops):
        #state:[[REG ID order],[CODE NO. order] ]
        state =[ [ random.randint(0,3) for i in range(Rnum)],[ random.randint(0,1) for i in range(Ilen)] ] 
        ca = CodeAnnealer(state,core_instlist_gen)
        ca.steps = 10000 * (Ilen//20)
        ca.Tmax = 100
        ca.Tmin = 1e-3
        res_state, _ = ca.anneal()
        res_instlist = getRefined(core_instlist_gen,res_state)
        res_li.append( [res_state, res_instlist, len(res_instlist) ] )
        #print (len(res_instlist), res_state)
    res_li.sort(key=lambda v:v[2])

    # generated_instlist, found_state
    return res_li[0][1], res_li[0][0]

#######################################################
def getRefined(instlist_orig,state, DEBUG=False):
    rsmn_li, t_li = state

    # utilitiy function!
    isImm = lambda reg: reg.split('#')[0] in ["VF","VI","C"]
    regtype = lambda reg: get_reg_name(reg)

    regtype_MNTRS = lambda reg: {"I":"M","J":"R","M":"M","N":"N","T":"T","R":"R","S":"S","F":"S","A":"A","U":"U","Z":"Z","O":"O"}[regtype(reg)]
    def check_RWconflict( reads, writes ):
        reads_wo_imm = [ reg for reg in reads if not isImm(reg) ]
        rtps = [ regtype_MNTRS(reg) for reg in reads_wo_imm ]
        wtps = [ regtype_MNTRS(reg) for reg in writes ]
        readM  = any([ rtp == 'M' for rtp in rtps ])
        writeM = any([ wtp == 'M' for wtp in wtps ])
        readN  = any([ rtp == 'N' for rtp in rtps ])
        writeN = any([ wtp == 'N' for wtp in wtps ])
        conflict_M = readM and writeM
        conflict_N = readN and writeN
        res = conflict_M or conflict_N
        return res

    # set state to instlist
    # mend 何もしない＆マスクの終了
    # mask 何もしない&マスクの開始
    # cmp 浮動小数点のsub。dvadd A -B してA-B>=0ならTrue。他のcmpは要対応

    num_unit = 2 # pe, alu
    unit_selector = {
        "add":0,
        "sub":0,
        "mul":0,
        "mulu":0,
        "muld":0,
        "fma":0,
        "fmau":0,
        "fmad":0,
        "fms":0,
        "fmsu":0,
        "fmsd":0,
        "fmi":0,
        "fmiu":0,
        "fmid":0,
        "set":1,
        "rsqrt":1,
        "max":1,
        "min":1,
        "nop:":0,
        "cmp":0,
        "lnot":1,
        "land":1,
        #"lor":1,
        "and":1,
        "xor":1,
        "dec":1,
        "rotate":1,
        "imm":1,
        "mask":2,
    }    
    instlist = []
    rotate_instlist = []
    maskdepth = 0
    empty_inst = ["",[],[], 0]
    for lno, line in enumerate(instlist_orig):
        op, reads, writes = line
        op_code, op_type = op.split(":")
        if op_code in ["to_f32","to_f32vec","to_f64","to_f64vec"]:
            print(f"error: {op_code} : type conversion is not supported for now!")
            assert False
        
        reads = [ "RSMN"[rsmn_li[int(get_reg_id(reg))]] + reg[1:] if regtype(reg) =="X" else reg for reg in reads]
        writes = [ "RSMN"[rsmn_li[int(get_reg_id(reg))]] + reg[1:] if regtype(reg)=="X" else reg for reg in writes]
        
        if op_code == "mend": # this is not real operator.
            maskdepth -= 1
            continue 
            
        if op_code == "cmp":
            if op_type[0] == "F":
                op_code = "sub"
                op = op_code + ":" + op_type
            else:
                op_code = "isub"
                op = op_code + ":" + op_type

        if op_code == "mask":
            maskdepth += 1
            continue

        if op_code == "rotate":
            if regtype_MNTRS(reads[0]) == regtype_MNTRS(reads[1]):
                return [ [] for i in range(100000) ]
            rotate_instlist.append([op,reads,writes,0])
            continue
        unit = unit_selector[op_code]
        # check if can be overlapped
        instbuf=[]
        while (len(instlist)>0):
            last_insts = instlist[-1]
            canOverlapped = True
            for r in reads + ([f"O.{maskdepth}:BOOL"] if maskdepth>0 else []):
                read = r.split(':')[0]
                for inst in last_insts:
                    if any( [ read == w.split(':')[0] for w in inst[2] ] ):
                        canOverlapped = False
            if canOverlapped:
                instbuf.append(instlist.pop())
            else:
                break
        isOverlapped = False
        while len(instbuf) > 0:
            if instbuf[-1][unit] == empty_inst:
                canOverlapped = True
                for inst in instbuf[-1]:
                    if op_code == "imm":
                        if any( [reg[0] in 'MI' for reg in inst[1]+inst[2]]):
                            canOverlapped = False
                if (not isOverlapped) and canOverlapped: 
                    instbuf[-1][unit] = [op,reads,writes,maskdepth]
                    isOverlapped = True
            instlist.append(instbuf.pop())
        if not isOverlapped:    
            instlist.append([ [op,reads,writes,maskdepth]  if i == unit else empty_inst for i in range(num_unit) ])

    if DEBUG:
        print_maskedinstlist( instlist)
        print ("-----refinement start----")
            
    #refine solver     
    
    iwc = {} # data race wait counter
    irc = {} # port conflict counter
    IWDIC = {"R":2,"S":2,"M":3,"N":3,"T":2,"A":0,"U":0,"Z":0,"O":0}
    IRDIC = {"R":0,"S":0,"M":3,"N":3,"T":0,"A":0,"U":0,"Z":0,"O":0}
    treg_cont = "" # current content of T-register


    remain_lines = [ it for it in instlist ][::-1]
    ret_instlist = []
    ret_len = 0
    #print("init_remain_lines",remain_lines)

    last_writes = [] # keep last wrote reg for forwarding

    SUCCESS = True
    save_reg_count = 0
    for i in range(10000):
    
        if len(remain_lines)==0:
            break
        for i in range(10000):
            if ret_len >= len(ret_instlist):
                break
            ret_len += 1
            iwc = { k:v-1 if v>0 else 0 for k,v in iwc.items() if v>0 } ## countdown of wait counter
            print (irc) if DEBUG else None
            irc = { k:v-1 if v>0 else 0 for k,v in irc.items() if v>0 } ## countdown of wait counter            
    
        line = remain_lines.pop()
        #print(line)
        reads = []
        writes = []
        optypes = []
        md = []
        for inst in line:
            op, r, w, md_ = inst
            reads += r
            writes += w
            optypes.append(op.split(":")[1]) if op != "" else None
            md.append(md_)

        #check conflict of reading
        def getAvailableReg(excludes):
            for regtype in ["R","S","M","N"]:
                if regtype in excludes:
                    pass
                else:
                    ret = regtype
                    return ret
            # can't reach
            raise AssertionError()
       
        # check noforward availability
        #can_noforward = True
        #if len(remain_lines) > 0:
        #    # check if already has noforward or nowrite
        #    has_noforward = False
        #    for inst in line:
        #        if inst[0] == "noforward:":
        #            has_noforward = True
        #    if not has_noforward:
        #        next_line = remain_lines.pop()
        #        for next_op, next_reads, next_writes, next_md in next_line:
        #            for w in writes:
        #                if w in next_reads:
        #                    can_noforward = False
        #        remain_lines.append(next_line)
        #        if can_noforward:
        #            line.append(["noforward:",[],[],md])

        #check if forwarding can apply
        if len(ret_instlist) > 0:
            forward_reads = reads
            for i in range(2): #assuming peinst and aluinst only
                prev_op, _, prev_writes, prev_md = ret_instlist[-1][i]
                if len(prev_op.split(":")) != 2:
                    continue
                pop, _ = prev_op.split(":")
                for pw in prev_writes:
                    if pw[0] == "O":
                        continue
                    if unit_selector[pop] == 0:
                        forward_reads = [ f"U.:{get_reg_type(reg)}" if reg == pw and all([m == prev_md or line[i][0] == "" for i,m in enumerate(md)]) else reg for reg in forward_reads ]
                    elif unit_selector[pop] == 1:
                        forward_reads = [ f"A.:{get_reg_type(reg)}" if reg == pw and all([m == prev_md or line[i][0] == "" for i,m in enumerate(md)]) else reg for reg in forward_reads ]
                    else:
                        print(f"unsupported operator {pop}")
                        assert False

            if reads != forward_reads and len(forward_reads)>0:
                # need to check conflicts here because forwarding can break if nop is inserted. Or, nop should always with noforward???
                port_conflict = False
                for reg in forward_reads:
                    if isImm(reg):
                        continue
                    reg_MNTRS = regtype_MNTRS(reg)
                    if irc.setdefault(reg_MNTRS,0) > 0:
                        port_conflict = True
                data_race = any([ iwc.setdefault(reg,0) > 0 for reg in forward_reads ]) 
                if not (port_conflict or data_race):
                    # get written vars in previous inst
                    prev_writes = []
                    for inst in ret_instlist[-1]:
                        prev_writes += inst[2]
                    # check nowrite availability
                    for w in prev_writes:
                        #print(w)
                        nowrite = True if ((get_reg_name(w) not in ["F","J","O"]) and (w not in forward_reads)) else False
                        if not nowrite:
                            continue
                        for l in reversed(remain_lines):
                            #print("checking",l)
                            finish = False
                            for op_,reads_,writes_,md_ in l:
                                if md_ > 0:
                                    reads.append(f'O.{md_}')
                                if len(reads_) == 0:
                                    continue
                                nowrite = nowrite and all([ reg.split(':')[0] != w.split(':')[0] for reg in reads_ ])
                                for w_ in writes_:
                                    #print(prev_writes[0],reads_,writes_,line,nowrite)
                                    if w_ == w:
                                        if any([ md_ > m for m in md ]):
                                            nowrite = False
                                        finish = True
                            if finish:
                                break
                        if nowrite:
                            for i in range(len(ret_instlist[-1])):
                                ret_instlist[-1][i][2] = [f"Z.:{get_reg_type(w_)}" if w_ == w else w_ for w_ in ret_instlist[-1][i][2]] 
                            iwc[prev_writes[0]] = 0
                    for i in range(len(ret_instlist[-1])):
                        if "Z." in ret_instlist[-1][i][2] and len(list(set(ret_instlist[-1][i][2]))) > 1:
                            ret_instlist[-1][i][2] = [ reg for reg in ret_instlist[-1][i][2] if get_reg_name(reg) != "Z" ]
                            assert len(ret_instlist[-1][i][2]) > 0

                    forward_reads.reverse()
                    ret=[]
                    for i in range(2):
                        r = []
                        for j in range(len(line[i][1])):
                            r.append(forward_reads.pop())
                        ret.append([line[i][0],r,line[i][2],line[i][3]])

                    remain_lines.append(ret)
                    #print("Forwarding applied")
                    continue

        # check read conflict
        tmp = reads[::-1]
        conf_list = []
        save_list = []
        while len(tmp) > 0:
            regA = tmp.pop()
            if isImm(regA):
                continue
            rtpA = regtype_MNTRS(regA)
            if rtpA in ['A','U']:
                continue
            isConflist = False
            for regB in tmp:
                if isImm(regB):
                    continue
                rtpB = regtype_MNTRS(regB)
                if rtpB in ['A','U']:
                    continue
                if ((rtpA == rtpB) and (regA != regB)):
                    conf_list.append(regA) if regA not in conf_list else None
                    isConflict = True
                    break
            if not (isConflist or (regA in save_list) or (regA in conf_list)):
                    save_list.append(regA)
        if len(conf_list) > 0: # conflict occur
            if len(conf_list)+len(save_list) > 4: # MNRS
                #print(line)
                #print("conf",conf_list)
                #print("save",save_list)
                for i in range(1000):
                    ret_instlist.append([])
                return ret_instlist
            replace_list = []
            for reg in conf_list:
                rtp = regtype_MNTRS(reg)
                save_reg = getAvailableReg([ regtype_MNTRS(r) for r in save_list ]) + f".s{save_reg_count}:{get_reg_type(reg)}" 
                save_reg_count += 1
                replace_list.append([reg,save_reg])
                save_list.append(save_reg)
            #print(replace_list)
            for replacer in replace_list:
                reads = [ replacer[1] if replacer[0] == r else r for r in reads ]
            reads.reverse()
            inst = []
            r = [[],[]]
            for i in range(2): # assuming peinst and aluinst only
                for j in range(len(line[i][1])):
                     r[i].append(reads.pop())
                inst.append([line[i][0],r[i],line[i][2],line[i][3]])
            remain_lines.append(inst)
            while any([ regtype_MNTRS(reg) in ['A','U'] for reg in remain_lines[-1][0][1] + remain_lines[-1][1][1] if not isImm(reg)]):
                remain_lines.append(ret_instlist.pop())
            for replacer in replace_list:
                tp = replacer[0].split(':')[1]
                remain_lines.append([empty_inst,[f"set:{tp}",[replacer[0]],[replacer[1]],0]])
            continue

        # check write conflict
        tmp = writes[::-1]
        conf_list = []
        save_list = []
        while len(tmp) > 0:
            regA = tmp.pop()
            rtpA = regtype_MNTRS(regA)
            if rtpA in ['Z']:
                continue
            isConflist = False
            for regB in tmp:
                rtpB = regtype_MNTRS(regB)
                if ((rtpA == rtpB) and (regA != regB)):
                    conf_list.append(regA) if regA not in conf_list else None
                    isConflict = True
                    break
            if not (isConflist or (regA in save_list) or (regA in conf_list)):
                    save_list.append(regA)
        if len(conf_list) > 0: # conflict occur
            #print(line)
            #print(conf_list)
            #print(save_list)
            assert len(conf_list)+len(save_list) <= 4 # MNRS
            replace_list = []
            for reg in conf_list:
                rtp = regtype_MNTRS(reg)
                tp = get_reg_type(reg)
                save_reg = getAvailableReg([ regtype_MNTRS(r) for r in save_list ]) + f".s{save_reg_count}:{tp}"
                save_reg_count += 1
                replace_list.append([reg,save_reg])
                save_list.append(save_reg)
            #print("replacer",replace_list)
            for replacer in replace_list:
                writes = [ replacer[1] if replacer[0] == w else w for w in writes ]
            writes.reverse()
            replaced_inst = []
            w = [[],[]]
            for i in range(len(line)):
                for j in range(len(line[i][2])):
                     w[i].append(writes.pop())
                replaced_inst.append([line[i][0],line[i][1],w[i],line[i][3]])
            for replacer in replace_list:
                for insts in reversed(remain_lines):
                    finish = False
                    for inst in insts:
                        #print(inst[1],replacer)
                        inst[1] = [ replacer[1] if r == replacer[0] else r for r in inst[1] ]
                        if any([ w == replacer[0] for w in inst[2] ]):
                            finish = True
                    if finish:
                        break
            remain_lines.append(replaced_inst)
            #print(remain_lines)
            #assert False
            continue
        
        # if RW conflict (ex. read N and write M)
        replace_list = []
        if check_RWconflict(reads,writes): #T.レジスタが変更されてても大丈夫か？ <- checkしてない？ by subarutaro
            wtps = [ regtype_MNTRS(w) for w in writes ]
            for w in writes:
                reg_MNTRS = regtype_MNTRS(w)
                if reg_MNTRS in "RSZO":
                    continue
                tp = get_reg_type(w)
                savereg_type = getAvailableReg(wtps)
                savereg = savereg_type + f".s{save_reg_count}:{tp}"
                save_reg_count += 1
                replace_list.append([w,savereg])
        if len(replace_list) > 0:
            for replacer in replace_list:
                writes = [ replacer[1] if replacer[0] == w else w for w in writes ]
            inst = []
            for i in range(2): # assuming peinst and aluinst only
                w = [[],[]]
                for j in range(len(line[i][2])):
                     w[i].append(writes.pop())
                inst.append([line[i][0],line[i][1],w[i],line[i][3]])
            for replacer in replace_list:
                remain_lines.append([empty_inst,["set:F64",[replacer[0]],[replacer[1]],0]])
            remain_lines.append(inst)
            continue
            
        #check port conflict
        port_conflict = False
        for reg in reads:
            if isImm(reg):
                continue
            reg_MNTRS = regtype_MNTRS(reg)
            if irc.setdefault(reg_MNTRS,0) > 0:
                port_conflict = True
        if port_conflict: 
            ret_instlist.append([["noforward:",[],[],0],empty_inst])
            remain_lines.append(line)
            print ("port_conflict_nop") if DEBUG else None
            continue    

        #check all readable and if not apply delay
        can_read = all([ iwc.setdefault(reg,0)==0 for reg in reads ])
        if not can_read:
            ret_instlist.append([["noforward:",[],[],0],empty_inst])
            remain_lines.append(line)
            print ("can_read_nop") if DEBUG else None
            continue
            
        # the other case (No regulation mathced, leave as it is)
        ret_instlist.append(line)
        if len(writes)>0:
            for write in writes:
                write_MNTRS = regtype_MNTRS(write)
                iwc[write] = IWDIC[write_MNTRS] 
                if len(write.split(':')) > 2: # simd case
                    name, tp, simd = write.split(':')
                    non_simd = name + ':' + tp
                    iwc[non_simd] = IWDIC[write_MNTRS] 
                else:
                    iwc[write+":0"] = IWDIC[write_MNTRS]
                    iwc[write+":1"] = IWDIC[write_MNTRS]


                irc[write_MNTRS] = IRDIC[write_MNTRS]
        print ("@",irc) if DEBUG else None

    # overlap rotate inst
    # rotate uses both input/output grf0 which tends to result in bad register allocation.
    # so, rotate insts are overlapped after refinement.
    for rotate in reversed(rotate_instlist):
        op,reads,writes,md = rotate
        isOverlapped = False
        regtype1 = regtype_MNTRS(reads[1])
        if regtype1 not in ["M","N"]: # need to check port conflict to use LM/LN
            for insts in reversed(ret_instlist):
                rs = []
                ws = []
                for inst in insts:
                    rs += [ reg for reg in inst[1] ]
                    ws += [ reg for reg in inst[2] ]
                if reads[0] in rs:
                    break
                if regtype1 in ["M","N"]:
                    if any([ regtype1 == regtype_MNTRS(reg) for reg in rs+ws ]):
                        continue
                if insts[unit_selector["rotate"]] != empty_inst:
                    continue
                canRead  = all([regtype_MNTRS(r) not in [ regtype_MNTRS(reg) for reg in rs ] for r in reads])
                canWrite = all([regtype_MNTRS(w) not in [ regtype_MNTRS(reg) for reg in ws ] for w in writes])
                if canRead and canWrite:
                    insts[unit_selector["rotate"]] = rotate
                    isOverlapped = True
                    break
        if not isOverlapped:
            ret_instlist.append([ rotate if i==unit_selector["rotate"] else empty_inst for i in range(num_unit) ])

    if SUCCESS:
        #print(ret_instlist)
        #sys.exit
        return ret_instlist
    else:
        if DEBUG:
            print ("Refinement failed. Ongoing result is following:")
            print_instlist(ret_instlist)
        return None

class CodeAnnealer(Annealer):
    def __init__(me, init_state, instlist_orig):
        super().__init__(init_state)
        me.instlist_orig = instlist_orig
        
    def move(me):
        rsmn_li, t_li = me.state
        reglen, instlen = len(rsmn_li), len(t_li)
        
        if random.randint(0,1)==0:
            #rsmn change
            ri = random.randint(0,reglen-1)
            rsmn_li[ ri ] = ( rsmn_li[ ri ] + 1 + random.randint(0,1) ) % 3
        else:
            #t change
            ri = random.randint(0,instlen-1)
            t_li[ ri ] = ( t_li[ ri ] + 1 ) % 2
        me.state = [rsmn_li, t_li]
    def energy(me):
        ret_instlist = getRefined(me.instlist_orig,me.state)
        rsmn_li, t_li = me.state
        mn_cost = sum([ 0.01 if regtype>0 else 0 for regtype in rsmn_li ])
        return len(ret_instlist) + mn_cost
