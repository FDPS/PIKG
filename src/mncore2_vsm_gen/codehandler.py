import subprocess as sb

def run_shell(cmd):
    try:
        out = sb.check_output(cmd,shell=True,stderr=sb.STDOUT)
        print (out.decode())
    except sb.CalledProcessError as exc:
        print("Status : FAIL", exc.returncode)
        print ("===")
        print( exc.output.decode() )
        
class Exe:
    def __init__(me,code):
        me.code = code
        sb.check_output('rm -rf ./tmp.*',shell=True,stderr=sb.STDOUT)
        with open("tmp.vsm","w") as f:
            f.write(code)

        try:
            sb.check_output('$ASM tmp.vsm > tmp.asm',shell=True,stderr=sb.STDOUT)
        except sb.CalledProcessError as exc:
            print("Status : FAIL", exc.returncode)
            print ("===")
            print( exc.output.decode() )

class Code:
    def __init__(me):
        me.ops= []
    #def __call__(op,mask=None,vlen=4):
        #me.ops.append( (op,mask,vlen) )
    def __getitem__(me,item):
        if isinstance(item, tuple):
            op,mask,vlen,sole = item
        else:
            op,mask,vlen,sole = item,"mask 0",4,False
        me.ops.append( [op.strip(" ").strip("\n").strip(" "),mask,vlen,sole] )
    def __add__(me,other):
        ret = Code()
        ret.ops = me.ops + other.ops
        return ret
        
def ToTextCode(code):
    text = ""
    for op,mask,vlen,sole in code.ops:
        #text+="vlen {0}\n".format(vlen)
        #text+="{0}\n".format(mask)
        text+="{0}\n".format(op) # assuming single write operand
    return text

def WriteTextCode(code,name='code'):
   # print ("generating code, target:", name)
    sb.check_output('rm -rf {0}.vsm'.format(name),shell=True,stderr=sb.STDOUT)
    with open("{0}.vsm".format(name),"w") as f:
        f.write(code)
    #run_shell('bash compile-addname.sh {0}'.format(name))
    #print ("original pack generated.")
    #with open("{0}.origpack".format(name)) as f:
    #    lines = f.readlines()
    #with open("{0}.pack".format(name),"w") as f:
    #    f.writelines([ line for line in lines if line[:2] in ['i ','m '] ])
    #print ("finished.")
    
def WriteTextCode_silent(code,name='code'):
    sb.run('rm -rf {0}.vsm'.format(name),shell=True,stderr=sb.STDOUT)
    with open("{0}.vsm".format(name),"w") as f:
        f.write(code)
    sb.check_output( 'bash compile-addname.sh {0}'.format(name) ,shell=True,stderr=sb.STDOUT)
    with open("{0}.origpack".format(name)) as f:
        lines = f.readlines()
    with open("{0}.pack".format(name),"w") as f:
        f.writelines([ line for line in lines if line[:2] in ['i ','m '] ])
    
