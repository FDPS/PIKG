def remove_type_conversion(instlist):
    dw_conversions = [ inst for inst in instlist if inst[0] in ["to_f32", "to_f32vec"] ]
    up_conversions = [ inst for inst in instlist if inst[0] in ["to_f64", "to_f64vec"] ]
    insts = [ inst for inst in instlist if inst[0] not in ["to_f32", "to_f32vec", "to_f64", "to_f64vec"] ]
    dw_convertor = {}
    for op, reads, writes in dw_conversions:
        assert len(reads) == 1 and len(writes) == 1
        read = reads[0]
        write = writes[0]
        if read[0] == 'M':
            dw_convertor[write] = read
        else:
            dw_convertor[read] = write
    up_convertor = { write[0]:read[0] for op,read,write in up_conversions }
    ret = []
    for inst in insts:
        op, reads, writes = inst
        if len(writes)>0:
            reads = [ up_convertor[read]  if read in up_convertor.keys() else read for read in reads ]
            writes = [ up_convertor[write]  if write in up_convertor.keys() else write for write in writes ]
            reads = [ dw_convertor[read]  if read in dw_convertor.keys() else read for read in reads ]
            writes = [ dw_convertor[write]  if write in dw_convertor.keys() else write for write in writes ]
            if any([ "64" in r.split(':')[1] for r in reads ]):
                op = op + ":" + writes[0].split(":")[1].replace("32","64")
            else:
                op = op + ":" + writes[0].split(":")[1]
        ret.append((op, reads,writes))
    return ret
