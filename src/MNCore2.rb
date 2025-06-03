require "json"

class Kernelprogram
  def generate_json_decllist(h=$varhash)
    decllist = []
    h.each{ |decl|
      if decl[1][0] =~ /(EPI|EPJ|FORCE)/
        name = decl[0]
        iotype, type, fdpsname, attr = decl[1]
        d = "decl_"
        d += "local_" if attr == "local"
        d += iotype.downcase
        fdpsname = iotype + "." + fdpsname if name =~ /(EPI_|EPJ_|FORCE_)/
        name = name.gsub(/EPI_/,'EPI.')
        name = name.gsub(/EPJ_/,'EPJ.')
        name = name.gsub(/FORCE_/,'FORCE.')
        decllist.append([d,[["TYPEIS",type], ["DECL_VAR",name], ["MEMBER_VAR",fdpsname]]])
      end
    }
    save_list = { "decllist" => decllist }
    File.write('./decllist.json',JSON.dump(save_list))
  end

  def generate_optimized_mncl_kernel_for_mncore2(filename,conversion_type,h = $varhash)
    #abort "error: --class-file option must be specified for conversion_type MNCore2" if $epi_file == nil

    generate_json_decllist(h)
    gen_packs_options = ""
    if $max_subwalk_size > 0
      gen_packs_options += " --max-subwalk-size #{$max_subwalk_size}"
    end
    ret = `mkdir -p packs && python3 #{File.dirname(__FILE__)}/mncore2_vsm_gen/gen_codehead.py #{filename} && python3 #{File.dirname(__FILE__)}/mncore2_vsm_gen/gen_packs.py #{gen_packs_options}`
    abort "error: gen vsm failed" if ret != "finished\n"

    regs_dic=[]
    decl_list=[]
    epi_json_list=[]
    epj_json_list=[]
    force_json_list=[]
    json_file="./codehead.json"
    File.open(json_file) { |file|
      hash = JSON.load(file)
      regs_dic = hash["regs_dic"]
      decl_list= hash["decllist"]
    }
    epi_json_list = regs_dic["I"]
    epj_json_list = regs_dic["J"]
    force_json_list = regs_dic["F"]
    epi_list = []
    epj_list = []
    force_list = []
    epi_json_list.each{ |v|
      v = v.gsub(/EPI./,"EPI_")
      epi_list.push(v.split(':')[0].split('@')[0])
    }
    epi_list.uniq!
    epj_json_list.each{ |v|
      v = v.gsub(/EPJ./,"EPJ_")
      epj_list.push(v.split(':')[0].split('@')[0])
    }
    epj_list.uniq!
    force_json_list.each{ |v|
      v = v.gsub(/FORCE./,"FORCE_")
      force_list.push(v.split(':')[0].split('@')[0])
    }
    force_list.uniq!
     
    fvars = []
    $max_element_size = 16
    $min_element_size = 64
    max_j_element_size = 16
    min_j_element_size = 64
    io_total_size = [0,0,0]
    io_element_type = [[],[],[]]
    ['I','J','F'].each_with_index{ |io,i|
      regs_dic[io].each{ |l|
        var, type = l.split(':')
        iotype, v = var.split('.')
        fdpsname = v.split('@')[0]
        v = iotype + '_' + fdpsname
        $max_element_size = [$max_element_size,get_single_data_size(type)].max
        $min_element_size = [$min_element_size,get_single_data_size(type)].min
        max_j_element_size = [max_j_element_size, get_single_data_size(type)].max
        min_j_element_size = [min_j_element_size, get_single_data_size(type)].min
        io_total_size[i] += get_data_size(type)
        io_element_type[i].push(type)
      
        if h[v] == nil
          abort "error: pikg source of pikg and mncore gencode are different!"
        end
        fvars.push(v) if not fvars.include?(v)
      }
    }
    $stderr.puts("worning: mixed precision in EPJ is not tested for DRAM indirect implementation") if max_j_element_size != min_j_element_size
    $max_inst_prec = 16
    $min_inst_prec = 64
    regs_dic.each{ |dic|
      dic[1].each{ |reg|
        tp = reg.split(":")[1]
        if tp == "BOOL"
            next
        end
        size = get_single_data_size(tp)
        $max_inst_prec = [$max_inst_prec,size].max
        $min_inst_prec = [$min_inst_prec,size].min
      }
    }
    is_multi_prec = ($max_inst_prec != $min_inst_prec)
    #fvars = generate_force_related_map(@statements)

    nepjsh = 1000
    
    code = String.new

    code += "#include \"#{$epi_file}\"\n" if $epi_file != nil
    code += "#include \"#{$epj_file}\"\n" if $epj_file != $epi_file
    code += "#include \"#{$force_file}\"\n" if ($force_file != $epi_file) && ($force_file != $epj_file)
    code += "#include \"pikg_mncore2_util.hpp\"\n"
    code += "#include \"pikg_vector.hpp\"\n"
    code += "#include \"pikg_time_profiler.hpp\"\n"
    code += "PIKG::TimeProfiler prof;\n"
    code += "#include<mncl/host/cl/cl.h>\n"
    code += "#include <map>\n"
    code += "#include <algorithm>\n"
    code += "#ifdef PIKG_ENABLE_PERFETTO_TRACE\n"
    code += "#include <codegen/base/perfetto.h>\n"
    code += "#include <codegen/base/perfetto_trace.h>\n"
    code += "#endif\n"

    code += "#define PIKG_MN_CHECK(STATUS) if(STATUS != CL_SUCCESS){ std::cout << \"error: at line\" << __LINE__ << \", PIKG_MN_CHECK_FAILURE\"  << std::endl; }\n"

    code += "#ifdef _OPENMP\n"
    code += "#include <omp.h>\n"
    code += "#else\n"
    code += "int omp_get_thread_num(){ return 0; }\n"
    code += "int omp_get_num_threads(){ return 1; }\n"
    code += "#endif\n"

    code += "cl_kernel load_kernel( cl_context& context, std::string filename, std::string kernelname ){\n"
    code += "  static std::map<std::string, cl_kernel> kernelcache;\n"
    code += "  if ( kernelcache.count(filename) == 1 ){\n"
    code += "    return kernelcache[filename];\n"
    code += "  }\n"
    code += "  cl_int status;\n"
    code += "  cl_kernel kernel = load_kernel_from_file(context,filename,kernelname,status); PIKG_MN_CHECK(status);\n"
    code += "  kernelcache[filename] = kernel;\n"
    code += "  return kernel;\n"
    code += "}\n"

    code += "//CONST VALUES\n"
    code += "constexpr int NGROUP = 4;\n"
    code += "constexpr int NL2BpGROUP = 2;\n"
    code += "constexpr int NL1BpL2B = 8;\n"
    code += "constexpr int NMABpL1B = 16;\n"
    code += "constexpr int NPEpMAB = 4;\n"

    code += "constexpr int NL2B = NGROUP * NL2BpGROUP; // 8 \n"
    code += "constexpr int NL1B = NL2B * NL1BpL2B;     // 64\n"
    code += "constexpr int NMAB = NL1B * NMABpL1B;     // 1024\n"
    code += "constexpr int NPE  = NMAB * NPEpMAB;      // 4096\n"
    code += "constexpr int NPEpL1B = NMABpL1B * NPEpMAB; // 64\n"

    #simd_width = 64/$min_element_size
    simd_width = 64/$min_inst_prec
    ni_per_l1b = ((((4096 - 256) / ([io_total_size[0],io_total_size[2]].max / 32) / 64)*64)/((simd_width)*64))*(simd_width)*64
    if $max_subwalk_size > 0 and $max_subwalk_size <= ni_per_l1b
        ni_per_l1b = simd_width * (($max_subwalk_size / 64) * 64)
    end
    code += "constexpr int NIpL1B = #{ni_per_l1b};\n"
    nj_per_pe = [256 / (io_total_size[1]/32), (4*1024*1024/8)/(2*(io_total_size[1]/32)*4096)].min
    nj_per_pe = (nj_per_pe/(4*simd_width))*(4*simd_width) # to distribute 4LW per l1bmd
    code += "constexpr int NJpPE = #{nj_per_pe};\n"
    code += "constexpr int BSIZE = NL1B;\n"
    code += "constexpr int NIMAX = NIpL1B;\n"
    code += "constexpr int NJMAX = NJpPE * NPEpL1B;\n"
    code += "constexpr int LW_PER_DAR_ENTRY = 16;\n"
    code += "constexpr int JBLOCK_SIZE = 4*#{64/$min_element_size};\n"
    code += "constexpr int NJBLOCKMAX = NJMAX/(#{64/$min_element_size}*LW_PER_DAR_ENTRY);\n"
    code += "constexpr int DRAM_OFFSET_SIZE = 4*1024*1024/8;\n"
    code += "constexpr int DRAM_IN_OFFSET = 0;\n"
    code += "constexpr int DRAM_OUT_OFFSET = DRAM_IN_OFFSET + DRAM_OFFSET_SIZE;\n"
    code += "constexpr int DRAM_IN_J_OFFSET = DRAM_OUT_OFFSET +  DRAM_OFFSET_SIZE;\n"
    code += "constexpr int DRAM_KERNEL_OFFSET = DRAM_IN_J_OFFSET + 2*DRAM_OFFSET_SIZE;\n"
    code += "constexpr int MAX_KERNEL_LAUNCH = 16;\n"

    code += "constexpr int DRAM_INDIRECT_OFFSET_SIZE = 32*1024*1024;\n"
    code += "constexpr int DRAM_J_INDIRECT_OFFSET = DRAM_KERNEL_OFFSET*MAX_KERNEL_LAUNCH;\n"

    code += "constexpr int l2bm_chunk_size = 64;\n"

    code += "constexpr int N_WALK_LIMIT = 1024;\n"
    code += "constexpr int N_EPI_MAX = 1024*1024;\n"

    code += "int ID_EPJ_DUMMY; // ID_EPJ_DUMMY must be set when send_flag = true in Index mode \n"
    code += "int DEBUG_TARGET = 0;\n"

    ["EPI","EPJ","FORCE"].each{ |io|
      code += "class #{io}Buffer{\n"
      code += "  public:\n"
      fvars.each{ |v|
        iotype = h[v][0]
        if iotype == io
          type = h[v][1]
          fdpsname = h[v][2]
          if type =~ /vec/
            get_vector_elements(type).each{ |s|
              code += "  PIKG::#{get_single_element_type(type)}* #{fdpsname}#{s};\n"
            }
          else  
            code += "  PIKG::#{type}* #{fdpsname};\n"
          end
        end
      }
      code += "  PIKG::F32* buffer;\n"
      code += "};\n"
    }
    code += "static bool init_call = true;\n"

    code += "struct ParticleData{\n"
    code += "  EPIBuffer hEPI[MAX_KERNEL_LAUNCH];\n"
    code += "  EPJBuffer hEPJ[MAX_KERNEL_LAUNCH];\n"
    code += "  FORCEBuffer hFORCE[MAX_KERNEL_LAUNCH];\n"
    code += "  PIKG::S32* hIdEPJ[MAX_KERNEL_LAUNCH];\n"
    code += "  int NJ_CLUSTER_MAX[MAX_KERNEL_LAUNCH];\n"
    code += "  int EPJ_INC_SIZE[MAX_KERNEL_LAUNCH];\n"
    code += "  void init_allocate(){\n"
    code += "    for(int tag=0;tag<MAX_KERNEL_LAUNCH;tag++){\n"
    code += "      NJ_CLUSTER_MAX[tag] = 4;\n"
    ["EPI","EPJ","FORCE"].each_with_index{ |io,i|
      size="NIMAX*BSIZE"
      size="NJMAX*BSIZE*NJ_CLUSTER_MAX[tag]" if io == "EPJ"
      code += "      h#{io}[tag].buffer = new PIKG::F32[#{io_total_size[i]/32}*#{size}];\n"
      code += "      PIKG::U32 #{io}_offset = 0;\n"
    }
    nelem_epj = 0
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      if iotype =~ /(EPI|EPJ|FORCE)/
        size="NIMAX*BSIZE"
        size="NJMAX*BSIZE" if iotype == "EPJ"
        get_vector_elements(type).each{ |s|
          code += "      h#{iotype}[tag].#{fdpsname}#{s} = (PIKG::#{get_single_element_type(type)}*)&h#{iotype}[tag].buffer[#{iotype}_offset];\n"
          code += "      #{iotype}_offset += #{get_single_data_size(type)/32}*#{size};\n"
          nelem_epj += 1 if iotype == "EPJ"
        }
      end
    }
    code += "      EPJ_INC_SIZE[tag] = EPJ_offset;\n"
    code += "      hIdEPJ[tag] = new PIKG::S32[#{nelem_epj}*NJMAX*BSIZE*NJ_CLUSTER_MAX[tag]];\n"
    code += "    }\n"
    code += "  }\n"
    
    code += "  void reallocate_epj (int tag, int NEW_CLUSTER_SIZE){\n"
    code += "    NJ_CLUSTER_MAX[tag] = NEW_CLUSTER_SIZE;\n"
    code += "    if(hEPJ[tag].buffer != nullptr) delete[]  hEPJ[tag].buffer;\n"
    code += "    hEPJ[tag].buffer = new PIKG::F32[#{io_total_size[1]/32}*NJMAX*BSIZE*NJ_CLUSTER_MAX[tag]];\n"
    code += "    PIKG::U32 offset = 0;\n"
    fvars.each{ |v|
        iotype = h[v][0]
        if iotype == "EPJ"
          type = h[v][1]
          fdpsname = h[v][2]
          get_vector_elements(type).each{ |s|
            code += "    hEPJ[tag].#{fdpsname}#{s} = (PIKG::#{get_single_element_type(type)}*)&hEPJ[tag].buffer[offset];\n"
            code += "    offset += #{get_single_data_size(type)/32}*NJMAX*BSIZE;\n"
          }
        end
    }
    code += "    EPJ_INC_SIZE[tag] = offset;\n"
    code += "    std::cout << \"Note: New Memory reallocation:\" << NJ_CLUSTER_MAX << std::endl;\n"
    code += "  }\n"
    code += "  void reset_epj_adr (int tag){\n"
    code += "    PIKG::U32 offset = 0;\n"
    fvars.each{ |v|
        iotype = h[v][0]
        if iotype == "EPJ"
          type = h[v][1]
          fdpsname = h[v][2]
          get_vector_elements(type).each{ |s|
            code += "    hEPJ[tag].#{fdpsname}#{s} = (PIKG::#{get_single_element_type(type)}*)&hEPJ[tag].buffer[offset];\n"
            code += "    offset += #{get_single_data_size(type)/32}*NJMAX*BSIZE;\n"
          }
        end
    }
    code += "    EPJ_INC_SIZE[tag] = offset;\n"
    code += "  }\n"
    code += "  ~ParticleData(){\n"
    code += "    for(int tag=0;tag<MAX_KERNEL_LAUNCH;tag++){\n"
    ["EPI","EPJ","FORCE"].each{ |iotype|
      code += "      if(h#{iotype}[tag].buffer != nullptr){\n"
      code += "        delete[]  h#{iotype}[tag].buffer;\n"
      code += "      }\n"
    }
    code += "    }\n"
    code += "  }\n"
    code += "  void reallocate_epj_for_indirect(int nepj){\n"
    code += "    if(hEPJ[0].buffer != nullptr) delete[] hEPJ[0].buffer;\n"
    code += "    hEPJ[0].buffer = new PIKG::F32[#{io_total_size[1]/32}*nepj];\n"
    code += "    int offset = 0;\n"
    fvars.each{ |v|
        iotype = h[v][0]
        if iotype == "EPJ"
          type = h[v][1]
          fdpsname = h[v][2]
          get_vector_elements(type).each{ |s|
            code += "    hEPJ[0].#{fdpsname}#{s} = (PIKG::#{get_single_element_type(type)}*)&hEPJ[0].buffer[offset];\n"
            code += "    offset += #{get_single_data_size(type)/32}*nepj;\n"
          }
        end
    }
    code += "  }\n"
    code += "  void reallocate_epj_index (int tag, int NEW_CLUSTER_SIZE){\n"
    code += "    NJ_CLUSTER_MAX[tag] = NEW_CLUSTER_SIZE;\n"
    code += "    if(hIdEPJ[tag] != nullptr) delete[] hIdEPJ[tag];\n"
    code += "    hIdEPJ[tag] = new PIKG::S32[NJMAX*BSIZE*NJ_CLUSTER_MAX[tag]];\n"
    code += "  }\n"
    code += "};\n"

    code += "static ParticleData pdatalist[16];\n"

    code += "typedef struct _subwalk{\n"
    code += "    int walk_id;\n"
    code += "    int ni_off;\n"
    code += "    int ni_size;\n"
    code += "    int nj_size;\n"
    code += "} Subwalk;\n"

    code += "void divide_walks (\n"
    code += "  std::vector<Subwalk> &swalks,\n"
    code += "  const PIKG::S32          n_walk,\n"
    code += "  const PIKG::S32          n_epi[],\n"
    code += "  const PIKG::S32          n_epj[] ){\n"

    code += "  prof.increaseLevel();\n"

    code += "  prof.start(\"MakeSubwalk\");\n"
    code += "  prof.increaseLevel();\n"
    code += "  for (int iw=0; iw< n_walk; iw++){\n"
    code += "    for (int ni_off=0; ni_off<n_epi[iw]; ni_off+=NIMAX){\n"

    code += "      // subwalk size determined\n"
    code += "      int ni_size = n_epi[iw]-ni_off;\n"
    code += "      ni_size = (ni_size>NIMAX)? NIMAX: ni_size;\n"
    code += "      int nj_size = n_epj[iw];\n"

    code += "      Subwalk subwalk{iw, ni_off, ni_size, nj_size"
    code += "};\n"

    code += "      swalks.push_back( subwalk );\n"
    code += "    }\n"
    code += "  }\n"
    code += "  prof.decreaseLevel();\n"
    code += "  prof.end(\"MakeSubwalk\");\n"
    code += "  prof.start(\"SortWalk\");\n"
    code += "  std::sort(swalks.begin(),swalks.end(), [](Subwalk const &a, Subwalk const &b){return a.nj_size > b.nj_size;} );\n"
    code += "  prof.end(\"SortWalk\");\n"
    code += "  prof.decreaseLevel();\n"
    code += "}\n"

    code += "int make_transfer_buffer(const int tag, const std::vector<Subwalk> &swalks, const #{$epi_name}* epi[], const #{$epj_name}* epj[], int sw_off, int thread_id ){\n"

    code += "  prof.increaseLevel();\n"

    code += "  //calculate size (NJ)\n"
    code += "  int nj_size_max = 0;\n"
    code += "  for (size_t sw=0; sw<BSIZE; sw++){\n"
    code += "    if (sw+sw_off >= swalks.size() ) break;\n"
    code += "    const Subwalk& walk = swalks[sw+sw_off];\n"
    code += "    nj_size_max = (walk.nj_size>nj_size_max)? walk.nj_size: nj_size_max;\n"
    code += "  }\n"
    code += "  int SIZEOF_NJ_CLUSTER = (nj_size_max + NJMAX - 1)/ NJMAX;\n"

    code += "  // reallocate memory if lack of memory size\n"
    code += "  if ( SIZEOF_NJ_CLUSTER > pdatalist[thread_id].NJ_CLUSTER_MAX[tag] ){\n"
    code += "   pdatalist[thread_id].reallocate_epj(tag,SIZEOF_NJ_CLUSTER * 2);\n"
    code += "  }\n"

    code += "  EPIBuffer&   hEPI   = pdatalist[thread_id].hEPI[tag];\n"
    code += "  EPJBuffer&   hEPJ   = pdatalist[thread_id].hEPJ[tag];\n"
    code += "  FORCEBuffer& hFORCE = pdatalist[thread_id].hFORCE[tag];\n"

    code += "  prof.start(\"CopyEPI\");\n"
    code += "  #pragma omp parallel for\n"
    code += "  for (int i=0; i<NIMAX*BSIZE; i++){\n"
    code += "    const int sw = i/NIMAX;\n"
    code += "    if ( sw+sw_off < swalks.size() ){\n"
    code += "      const auto& walk = swalks[sw+sw_off];\n"
    code += "      const int w = walk.walk_id;\n"
    code += "      const int ni_off = walk.ni_off;\n"
    code += "      const int ii = i%NIMAX;\n"
    if $min_element_size == 32
      code += "      const int index32 = (BSIZE*(ii/(2*l2bm_chunk_size)) + sw)*(2*l2bm_chunk_size) + ((ii/2)%l2bm_chunk_size)*2 + (1^(ii%2));\n"
    end
    if $max_element_size == 64
      code += "      const int index64 = (BSIZE*(ii/l2bm_chunk_size) + sw)*l2bm_chunk_size + ii%l2bm_chunk_size;\n"
    end
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      size = get_single_data_size(type)
      if iotype == "EPI"
        get_vector_elements(type).each{ |s|
          code += "      h#{iotype}.#{fdpsname}#{s}[index#{size}] = epi[w][ni_off + ii].#{fdpsname}"
          code += ".#{s}" if s != ""
          code += ";\n"
        }
      end
    }
    code += "    }\n"
    code += "  }\n"
    code += "  prof.end(\"CopyEPI\");\n"

    code += "  prof.start(\"CopyEPJ\");\n"
    code += "  #pragma omp parallel for\n"
    code += "  for (int j=0; j<NJMAX*BSIZE; j++){\n"
    code += "    const int sw = j/NJMAX;\n"
    code += "    if(sw+sw_off < swalks.size()){\n"
    code += "      const auto& walk = swalks[sw+sw_off];\n"
    code += "      const int w = walk.walk_id;\n"
    code += "      const int jj = j%NJMAX;\n"
    # EPJ is distributed by l1bmd. each pe has data in order of 64*j + $peid
    if $min_element_size == 32
        code += "      const int index32 = ((jj/(2*l2bm_chunk_size))*BSIZE + sw)*(2*l2bm_chunk_size) + (jj%(l2bm_chunk_size))*2 + (1^((jj/l2bm_chunk_size)%2));\n"
    end
    if $max_element_size == 64
      code += "      const int index64 = ((jj/l2bm_chunk_size)*BSIZE + sw)*l2bm_chunk_size + jj%l2bm_chunk_size;\n"
    end
    code += "      const #{$epj_name} epj_dummy;\n"
    code += "      if(jj < walk.nj_size){\n"
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      size = get_single_data_size(type)
      if iotype == "EPJ"
        get_vector_elements(type).each{ |s|
          code += "        h#{iotype}.#{fdpsname}#{s}[index#{size}] = epj[w][jj].#{fdpsname}"
          code += ".#{s}" if s != ""
          code += ";\n"
        }
      end
    }
    code += "      }else{\n"
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      size = get_single_data_size(type)
      if iotype == "EPJ"
        get_vector_elements(type).each{ |s|
          if s == ""
            code += "        h#{iotype}.#{fdpsname}#{s}[index#{size}] = epj_dummy.#{fdpsname};\n"
          else
            code += "        h#{iotype}.#{fdpsname}#{s}[index#{size}] = epj_dummy.#{fdpsname}.#{s};\n"
          end
        }
      end
    }
    code += "      }\n"
    code += "    }\n"
    code += "  }\n"
    code += "  prof.end(\"CopyEPJ\");\n"

    code += "  prof.decreaseLevel();\n"
    code += "  return SIZEOF_NJ_CLUSTER;\n"
    code += "}\n"

    code += "static #{$force_name} **store;\n"
    code += "static std::vector<Subwalk> swalks;\n"
    code += "void retrieve_transfer_buffer( int tag, std::vector<Subwalk> &swalks, int sw_off, int thread_id){\n"
    code += "  FORCEBuffer hFORCE = pdatalist[thread_id].hFORCE[tag];\n"
    code += "  int ni_max = 0;\n"
    code += "  const int n_sw = std::min(BSIZE,(int)swalks.size()-sw_off);\n" # use NI_MAX instead?
    code += "  for(int sw=0; sw<n_sw; sw++) ni_max = std::max(ni_max,swalks[sw+sw_off].ni_size);\n"
    code += "  #pragma omp parallel for\n"
    code += "  for (size_t ii=0; ii<ni_max*n_sw; ii++){\n"
    code += "    const int sw = ii / ni_max;\n"
    code += "    const int i = ii % ni_max;\n"
    code += "    const Subwalk& walk = swalks[sw+sw_off];\n"
    code += "    const int w = walk.walk_id;\n"
    code += "    const int ni_off = walk.ni_off;\n"
    code += "    if(i < walk.ni_size){\n"
    if $min_element_size == 32
      code += "      const int index32 = NIMAX*sw + (i/2)*2 + (1^(i%2));\n"
    end
    if $max_element_size == 64
      code += "      const int index64 = NIMAX*sw + i;\n"
    end
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      size = get_single_data_size(type)
      if iotype == "FORCE"
        get_vector_elements(type).each{ |s|
          code += "        store[w][ni_off + i].#{fdpsname}"
          code += ".#{s}" if s != ""
          code += " = h#{iotype}.#{fdpsname}#{s}[index#{size}];\n"
        }
      end
    }
    code += "    }\n"
    code += "  }\n"
    code += "}\n"

    code += "cl_context context;\n"
    code += "cl_platform_id platform;\n"
    code += "cl_device_id devices[8];\n"
    code += "cl_uint ndevice;\n"
    code += "cl_command_queue queue;\n"
    code += "cl_program program;\n"
    code += "cl_kernel kernel;\n"
    nloop_max = 1024
    code += "constexpr int NLOOP_MAX = #{nloop_max};\n"
    code += "cl_event head_kernel_event[MAX_KERNEL_LAUNCH];\n"
    code += "cl_event kernel_events[MAX_KERNEL_LAUNCH][NLOOP_MAX];\n"
    code += "cl_event tail_kernel_event[MAX_KERNEL_LAUNCH];\n"
    code += "cl_event write_events[MAX_KERNEL_LAUNCH][1+NLOOP_MAX];\n"
    code += "cl_int status;\n"
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      if iotype =~ /EPJ/
        get_vector_elements(type).each{ |s|
          if iotype == "EPJ"
            code += "cl_mem #{iotype}#{fdpsname}#{s}BufferB0[MAX_KERNEL_LAUNCH];\n"
            code += "cl_mem #{iotype}#{fdpsname}#{s}BufferB1[MAX_KERNEL_LAUNCH];\n"
          else
            code += "cl_mem #{iotype}#{fdpsname}#{s}Buffer[MAX_KERNEL_LAUNCH];\n"
          end
        }
      end
    }
    code += "cl_mem dEPI[MAX_KERNEL_LAUNCH];\n"
    code += "cl_mem dEPJB0[MAX_KERNEL_LAUNCH],dEPJB1[MAX_KERNEL_LAUNCH];\n"
    code += "cl_mem dFORCE[MAX_KERNEL_LAUNCH];\n"
    code += "cl_mem DummyBuffer;\n"
    
    code += "int Initialize_this_header(){\n"
    code += "#ifdef PIKG_ENABLE_PERFETTO_TRACE\n"
    code += "  grape_pfn::base::StartTrace(\"perfetto_outtrace.pb\");\n"
    code += "#endif\n"
    code += "  std::cout << \"Initialized.\" << std::endl;\n"
    code += "  constexpr int device_to_use = 0;\n"
    code += "  status = clGetPlatformIDs(1, &platform, nullptr); PIKG_MN_CHECK(status);\n"
    code += "#ifdef USE_MNCORE2_EMULATOR // currently simulator is not supported\n"
    code += "  constexpr cl_device_type device_type = CL_DEVICE_TYPE_EMU_MNCORE2;\n"
    code += "#else\n"
    code += "  constexpr cl_device_type device_type = CL_DEVICE_TYPE_ACCELERATOR;\n"
    code += "#endif\n"
    code += "  status = clGetDeviceIDs(platform,device_type,0,nullptr,&ndevice); PIKG_MN_CHECK(status);\n"
    code += "  assert(ndevice>0 && ndevice<=8);\n"
    code += "  assert(device_to_use<ndevice);\n"
    code += "  status = clGetDeviceIDs(platform,device_type,ndevice,devices,nullptr); PIKG_MN_CHECK(status);\n"
    code += "  context = clCreateContext(nullptr, 1, &devices[device_to_use], nullptr, nullptr, &status); PIKG_MN_CHECK(status);\n"
    code += "  queue = clCreateCommandQueue(context, devices[device_to_use], 0, &status); PIKG_MN_CHECK(status);\n"
    code += "  status = clGetPlatformIDs(1, &platform, nullptr);PIKG_MN_CHECK(status);\n"
    code += "  char buffer[1024];\n"
    code += "  status = clGetPlatformInfo(platform, CL_PLATFORM_NAME, sizeof(buffer) - 1, buffer, NULL); PIKG_MN_CHECK(status);\n"
    code += "  std::cout << \"Platform Name: \" << buffer << std::endl;\n"

    code += "  status = clGetPlatformInfo(platform, CL_PLATFORM_VERSION, sizeof(buffer) - 1, buffer, NULL); PIKG_MN_CHECK(status);\n"
    code += "  std::cout << \"Platform Version: \" << buffer << std::endl;\n"

    code += "  status = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, sizeof(buffer) - 1, buffer, NULL);PIKG_MN_CHECK(status);\n"
    code += "  std::cout << \"Platform Vendor: \" << buffer << std::endl;\n"
    code += "  for(int i=0;i<ndevice;i++){\n"

    code += "    status = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, sizeof(buffer) - 1, buffer, NULL);PIKG_MN_CHECK(status);\n"
    code += "    std::cout << \"Device\"<<i<<\" Name: \" << \" \" << buffer << std::endl;\n"
    code += "  }\n"
    code += "  std::cout << \"Device\" << device_to_use << \" is used\" << std::endl;\n"

    #code += "  #pragma omp parallel\n"
    code += "  {\n"
    code += "    int thread_id = omp_get_thread_num();\n"
    code += "    pdatalist[thread_id].init_allocate();\n"
    code += "    std::cout << \"INIT ALLOCATION FOR tid= \" << thread_id << std::endl;\n"
    code += "  }\n"
    code += "  DummyBuffer = clCreateBufferWithAttributes(context, CL_MEM_READ_ONLY, 0, 0, sizeof(PIKG::F64)*64, nullptr, &status); PIKG_MN_CHECK(status);\n"
    code += "  return 0;\n"
    code += "}\n"

    code += "void waitForEventsByDummyRead(cl_command_queue& queue, int num_event, cl_event *event){\n"
    code += "  double buf[64];\n"
    code += "  status = clEnqueueReadBuffer(queue, DummyBuffer, CL_TRUE, 0, sizeof(PIKG::F64)*16, buf, num_event, event, nullptr); PIKG_MN_CHECK(status);\n"
    code += "}\n"

    code += "void DEVICE_RUN_DIVJ(int tag, int Nloop, int NI, int NJpPE, int thread_id, const std::vector<Subwalk> &swalks, const #{$epj_name}* epj[], int sw_off){ // Nloop (# of J accumeration)\n"

    code += "  EPIBuffer&   hEPI   = pdatalist[thread_id].hEPI[tag];\n"
    code += "  FORCEBuffer& hFORCE = pdatalist[thread_id].hFORCE[tag];\n"
    code += "  EPJBuffer&   hEPJ   = pdatalist[thread_id].hEPJ[tag];\n"

    code += "  prof.start(\"LoadKernels\");\n"
    code += "  cl_kernel Hkernel  = load_kernel(context,\"packs/nbY1_H\" + std::to_string(tag) + \".vsm\",\"Hkernel\" + std::to_string(tag));\n"
    code += "  cl_kernel FFkernel = load_kernel(context,\"packs/nbY1_MS\" + std::to_string(tag) + \".vsm\",\"FFkernel\" + std::to_string(tag));\n"
    code += "  cl_kernel RRkernel = load_kernel(context,\"packs/nbY1_MSR\" + std::to_string(tag) + \".vsm\",\"RRkernel\" + std::to_string(tag));\n"
    code += "  cl_kernel FTkernel  = load_kernel(context,\"packs/nbY1_MT\" + std::to_string(tag) + \".vsm\",\"FTkernel\" + std::to_string(tag));\n"
    code += "  cl_kernel RTkernel  = load_kernel(context,\"packs/nbY1_MRT\" + std::to_string(tag) + \".vsm\",\"RTkernel\" + std::to_string(tag));\n"
    code += "  prof.end(\"LoadKernels\");\n"

    code += "  // NI SET\n"
    code += "  prof.start(\"WriteEPI\");\n"
    code += "  status = clEnqueueWriteBuffer(queue, dEPI[tag], CL_FALSE, 0, #{io_total_size[0]/32}*sizeof(PIKG::F32)*NIMAX*BSIZE, hEPI.buffer, 0, nullptr, &write_events[tag][0]); PIKG_MN_CHECK(status);\n"
    code += "  prof.end(\"WriteEPI\");\n"
    code += "  // NJ SET\n"

    code += "  prof.start(\"WriteEPJFirst\");\n"
    code += "  status = clEnqueueWriteBuffer(queue, dEPJB0[tag], CL_FALSE, 0, #{io_total_size[1]/32}*sizeof(PIKG::F32)*NJMAX*BSIZE, hEPJ.buffer, 0, nullptr, &write_events[tag][1]); PIKG_MN_CHECK(status);\n"
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
       fdpsname = h[v][2]
       stype = get_single_element_type(type)
      if iotype == "EPJ"
        offset = "NJMAX*BSIZE*(loop+1)"
        size = "sizeof(PIKG::#{stype})*NJMAX*BSIZE"
        get_vector_elements(type).each{ |s|
          code += "      h#{iotype}.#{fdpsname}#{s} = (PIKG::#{stype}*)((PIKG::F32*)h#{iotype}.#{fdpsname}#{s} + pdatalist[thread_id].EPJ_INC_SIZE[tag]);\n"
        }
      end
    }
    code += "  prof.end(\"WriteEPJFirst\");\n"
    code += "  prof.start(\"EnqueueHeadKernel\");\n"
    code += "  status = clEnqueueTask(queue, Hkernel, 2, write_events[tag], &head_kernel_event[tag]); PIKG_MN_CHECK(status);\n"
    code += "  prof.end(\"EnqueueHeadKernel\");\n"

    #code += "#ifdef PIKG_MNCORE2_DEBUG\n"
    #code += "  cl_kernel debug_kernel = load_kernel(context, \"packs/nbY1_DBGH.vsm\",\"DBGkernel\");\n" 
    #code += "  cl_event debug_event;\n"
    #code += "  status = clEnqueueTask(queue,Hkernel, 1, &head_kernel_event[tag], &debug_event); PIKG_MN_CHECK(status);\n"
    #code += "  status = clEnqueueReadBuffer(queue, dFORCE[0], CL_TRUE, 0, #{io_total_size[1]/32}*sizeof(PIKG::F32)*NIMAX*BSIZE, hEPJ.buffer, 1, &debug_event, nullptr); PIKG_MN_CHECK(status);\n"
    #code += "  for(int i=0;i<NJMAX*BSIZE;i++){\n"
    #code += "  }\n"
    #code += "#endif\n"
   
    code += "  assert(Nloop <= NLOOP_MAX);\n" 
    code += "  for ( int loop=0; loop< Nloop ; loop++ ){\n"

    code += "    // NJ SET\n"
    code += "    if (loop<Nloop-1){\n"
    code += "#ifdef PIKG_PROFILE_COPY_EP_TO_BUFFER\n"
    code += "      auto s = std::chrono::high_resolution_clock::now();\n"
    code += "#endif\n"
    code += "      prof.start(\"CopyEPJtoBuffer\");\n"
    code += "      const int offset = NJMAX*BSIZE*(loop+1);\n"
    code += "      #pragma omp parallel for\n"
    code += "      for (int j=0; j<NJMAX*BSIZE; j++){\n"
    code += "        const int sw = j/NJMAX;\n"
    code += "        if(sw+sw_off < swalks.size()){\n"
    code += "          const auto& walk = swalks[sw+sw_off];\n"
    code += "          const int w = walk.walk_id;\n"
    code += "          const int jj = j%NJMAX;\n"
    if $min_element_size == 32
      code += "          const int index32 = ((jj/(2*l2bm_chunk_size))*BSIZE + sw)*(2*l2bm_chunk_size) + (jj%(l2bm_chunk_size))*2 + (1^((jj/l2bm_chunk_size)%2));\n"
    end
    if $max_element_size == 64
      code += "          const int index64 = ((jj/l2bm_chunk_size)*BSIZE + sw)*l2bm_chunk_size + jj%l2bm_chunk_size;\n"
    end
    code += "          const int jjj = jj + (loop+1)*NJMAX;\n"
    code += "          const #{$epj_name} epj_dummy;\n"
    code += "          if(jjj < walk.nj_size){\n"
    epj_total_size = 0
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      size = get_single_data_size(type)
      if iotype == "EPJ"
        get_vector_elements(type).each{ |s|
          code += "            h#{iotype}.#{fdpsname}#{s}[index#{size}] = epj[w][jjj].#{fdpsname}"
          code += ".#{s}" if s != ""
          code += ";\n"
          epj_total_size += size/8
        }
      end
    }
    code += "          }else{\n"
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      size = get_single_data_size(type)
      if iotype == "EPJ"
        get_vector_elements(type).each{ |s|
          if s == ""
            code += "            h#{iotype}.#{fdpsname}#{s}[index#{size}] = epj_dummy.#{fdpsname};\n"
          else
            code += "            h#{iotype}.#{fdpsname}#{s}[index#{size}] = epj_dummy.#{fdpsname}.#{s};\n"
          end
        }
      end
    }
    code += "          }\n"
    code += "        }\n"
    code += "      }\n"
    code += "      prof.end(\"CopyEPJtoBuffer\");\n"
    code += "#ifdef PIKG_PROFILE_COPY_EP_TO_BUFFER\n"
    code += "      auto e = std::chrono::high_resolution_clock::now();\n"
    code += "      auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(e-s).count() * 1e-9; // in sec\n"
    code += "      std::cout << \"ElapsedTime(sec): \" << elapsed << \" Bandwidth(GB/s): \" << 2*NJMAX*BSIZE*#{epj_total_size}/elapsed * 1e-9 << std::endl;\n"
    code += "#endif\n"
    code += "      cl_event& wait_kernel_event = (loop==0) ? head_kernel_event[tag] : kernel_events[tag][loop-1];\n"
    code += "      prof.start(\"WriteEPJ\");\n"
    nelem = "#{io_total_size[1]/32}*NJMAX*BSIZE"
    code += "      cl_mem& dEPJ = (loop%2==0) ? dEPJB1[tag] : dEPJB0[tag];\n"
    code += "      status = clEnqueueWriteBuffer(queue, dEPJ, CL_FALSE, 0, sizeof(PIKG::F32)*#{nelem}, hEPJ.buffer+#{nelem}*(loop+1), 1, &wait_kernel_event, &write_events[tag][1 + (loop+1)]); PIKG_MN_CHECK(status);\n"
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
       fdpsname = h[v][2]
       stype = get_single_element_type(type)
      if iotype == "EPJ"
        offset = "NJMAX*BSIZE*(loop+1)"
        size = "sizeof(PIKG::#{stype})*NJMAX*BSIZE"
        get_vector_elements(type).each{ |s|
          code += "      h#{iotype}.#{fdpsname}#{s} = (PIKG::#{stype}*)((PIKG::F32*)h#{iotype}.#{fdpsname}#{s} + pdatalist[thread_id].EPJ_INC_SIZE[tag]);\n"
        }
      end
    }
    code += "      prof.end(\"WriteEPJ\");\n"
    code += "    }\n"
    code += "    prof.start(\"EnqueueMainKernel\");\n"
    code += "    if(loop==Nloop-1){\n"
    code += "      cl_kernel& kernel = (loop%2==0)? FTkernel : RTkernel;\n"
    code += "      if(loop == 0){\n"
    code += "        status = clEnqueueTask(queue, kernel, 1, &head_kernel_event[tag], &tail_kernel_event[tag]); PIKG_MN_CHECK(status);\n"
    code += "      }else{\n"
    code += "        status = clEnqueueTask(queue, kernel, 1, &kernel_events[tag][loop-1], &tail_kernel_event[tag]); PIKG_MN_CHECK(status);\n"
    code += "      }\n"
    code += "    }else{\n"
    code += "      cl_kernel& kernel = (loop%2==0)? FFkernel : RRkernel;\n"
    code += "      status = clEnqueueTask(queue, kernel, 1, &write_events[tag][1+(loop+1)], &kernel_events[tag][loop]); PIKG_MN_CHECK(status);\n"
    code += "    }\n"
    code += "    prof.end(\"EnqueueMainKernel\");\n"
    code += "  }\n"
    code += "  pdatalist[thread_id].reset_epj_adr(tag);\n"
    code += "}\n"

    code += "void RUN_WALKS(\n"
    code += "                             #{$force_name}     *store[],\n"
    code += "                             const PIKG::S32      n_walk,\n"
    code += "                             const #{$epi_name} *epi[],\n"
    code += "                             const PIKG::S32      n_epi[],\n"
    code += "                             const #{$epj_name} *epj[],\n"
    code += "                             const PIKG::S32      n_epj[]){\n"

    code += "    if(init_call){\n"
    code += "                Initialize_this_header();\n"
    code += "                init_call = false;\n"
    code += "    }\n"

    code += "    if ((NJMAX%64)!=0) { std::cout << \"NJ != 0 mod 64 Error!\" << std::endl; exit(0); }\n"
    code += "    prof.start(\"RunWalks\");\n"
    code += "    prof.increaseLevel();\n"

    code += "    swalks.clear();\n"
    code += "    prof.start(\"DevideWalks\");\n"
    code += "    divide_walks( swalks, n_walk, n_epi, n_epj);\n"
    code += "    prof.end(\"DevideWalks\");\n"

    code += "    for (size_t sw_off=0, tag=0; sw_off< swalks.size(); sw_off+=BSIZE, tag++){\n"
    code += "      assert(tag < MAX_KERNEL_LAUNCH);\n"
    code += "      int thread_id = omp_get_thread_num();\n"

    code += "      prof.start(\"MakeTransferBuffer\");\n"
    code += "      int SIZEOF_NJ_CLUSTER = make_transfer_buffer(tag,swalks, epi, epj, sw_off, thread_id);\n"
    code += "      prof.end(\"MakeTransferBuffer\");\n"
    code += "      prof.start(\"CreateMNCLBuffer\");\n"
    code += "      {\n"
    code += "        int epi_offset   = sizeof(PIKG::F64)*(tag*DRAM_KERNEL_OFFSET + DRAM_IN_OFFSET);\n"
    code += "        int epj_offset   = sizeof(PIKG::F64)*(tag*DRAM_KERNEL_OFFSET + DRAM_IN_J_OFFSET);\n"
    code += "        int force_offset = sizeof(PIKG::F64)*(tag*DRAM_KERNEL_OFFSET + DRAM_OUT_OFFSET);\n"
    #fvars.each{ |v|
    #  iotype   = h[v][0]
    #  type     = h[v][1]
    #  fdpsname = h[v][2]
    #  stype = get_single_element_type(type)
    #  if iotype =~ /EPJ/
    #    read_or_write = "CL_MEM_WRITE_ONLY"
    #    read_or_write = "CL_MEM_READ_ONLY" if iotype == "FORCE"
    #    size = "sizeof(PIKG::#{stype})*NIMAX*BSIZE"
    #    size = "sizeof(PIKG::#{stype})*NJMAX*SIZEOF_NJ_CLUSTER" if iotype == "EPJ"
    #    increment = size
    #    increment = "sizeof(PIKG::#{stype})*NJMAX*BSIZE" if iotype == "EPJ"
    #    suffixes = [""]
    #    suffixes = ["B0","B1"] if iotype == "EPJ"
    #    get_vector_elements(type).each{ |s|
    #      offset = "epi_offset" if iotype == "EPI"
    #      offset = "epj_offset" if iotype == "EPJ"
    #      offset = "force_offset" if iotype == "FORCE"
    #      suffixes.each{ |suffix|
    #        isx = ""
    #        isx = " + sizeof(PIKG::F64)*DRAM_OFFSET_SIZE" if suffix == "B1"
    #        code += "        #{iotype}#{fdpsname}#{s}Buffer#{suffix}[tag] = clCreateBufferWithAttributes(context, #{read_or_write}, 0, #{offset}#{isx}, #{size}, nullptr, &status); PIKG_MN_CHECK(status);\n"
    #      }
    #      code += "        #{offset} += #{increment};\n"
    #    }
    #  end
    #}
    code += "        dEPI[tag] = clCreateBufferWithAttributes(context, CL_MEM_READ_ONLY, 0, sizeof(PIKG::F64)*(tag*DRAM_KERNEL_OFFSET+DRAM_IN_OFFSET), #{io_total_size[0]/32}*sizeof(PIKG::F32)*NIMAX*BSIZE, nullptr, &status); PIKG_MN_CHECK(status);\n"
    code += "        dEPJB0[tag] = clCreateBufferWithAttributes(context, CL_MEM_READ_ONLY, 0, sizeof(PIKG::F64)*(tag*DRAM_KERNEL_OFFSET+DRAM_IN_J_OFFSET), #{io_total_size[1]/32}*sizeof(PIKG::F32)*NJMAX*BSIZE*SIZEOF_NJ_CLUSTER, nullptr, &status); PIKG_MN_CHECK(status);\n"
    code += "        dEPJB1[tag] = clCreateBufferWithAttributes(context, CL_MEM_READ_ONLY, 0, sizeof(PIKG::F64)*(tag*DRAM_KERNEL_OFFSET+DRAM_IN_J_OFFSET+DRAM_OFFSET_SIZE), #{io_total_size[1]/32}*sizeof(PIKG::F32)*NJMAX*BSIZE*SIZEOF_NJ_CLUSTER, nullptr, &status); PIKG_MN_CHECK(status);\n"
    code += "        dFORCE[tag] = clCreateBufferWithAttributes(context, CL_MEM_READ_ONLY, 0, sizeof(PIKG::F64)*(tag*DRAM_KERNEL_OFFSET+DRAM_OUT_OFFSET), #{io_total_size[2]/32}*sizeof(PIKG::F32)*NIMAX*BSIZE, nullptr, &status); PIKG_MN_CHECK(status);\n"
    code += "      }\n"
    code += "      prof.end(\"CreateMNCLBuffer\");\n"

    code += "      prof.start(\"DeviceRunDivJ\");\n"
    code += "      prof.increaseLevel();\n"
    code += "      DEVICE_RUN_DIVJ( tag, SIZEOF_NJ_CLUSTER, NIMAX, NJpPE, thread_id, swalks, epj, sw_off);\n"
    code += "      prof.decreaseLevel();\n"
    code += "      prof.end(\"DeviceRunDivJ\");\n"

    code += "    }\n"
    code += "    prof.decreaseLevel();\n"
    code += "    prof.end(\"RunWalks\");\n"
    code += "}\n"

    code += "PIKG::S32 Dispatch#{$kernel_name}(\n"
    code += "                             const PIKG::S32          tag,\n"
    code += "                             const PIKG::S32          n_walk,\n"
    code += "                             const #{$epi_name}** __restrict__   epi,\n"
    code += "                             const PIKG::S32*    __restrict__   n_epi,\n"
    code += "                             const #{$epj_name}** __restrict__   epj,\n"
    code += "                             const PIKG::S32*    __restrict__   n_epj){\n"
    code += "    prof.start(\"Dispatch#{$kernel_name}\");\n"
    code += "    prof.increaseLevel();\n"
    code += "    prof.start(\"PrepareStoreBuffer\");\n"
    code += "    assert(n_walk<=N_WALK_LIMIT);\n"
    code += "    assert(n_epi[0]>0 && n_epj[0]>0);\n"
    code += "    static bool first = true;\n"
    code += "    if(first){\n"
    code += "      store = new #{$force_name}*[N_WALK_LIMIT];\n"
    code += "      for (int iw=0; iw<N_WALK_LIMIT; iw++){\n"
    code += "        store[iw] = new #{$force_name}[N_EPI_MAX];\n"
    code += "      }\n"
    code += "      first = false;\n"
    code += "    }\n"
    code += "    prof.end(\"PrepareStoreBuffer\");\n"

    code += "    RUN_WALKS( store, n_walk, epi, n_epi, epj, n_epj);\n"
    code += "    prof.decreaseLevel();\n"
    code += "    prof.end(\"Dispatch#{$kernel_name}\");\n"
    code += "    return 0;\n"
    code += "}\n"

    code += "PIKG::S32 Retrieve#{$kernel_name}(const PIKG::S32 tag,\n"
    code += "                       const PIKG::S32 n_walk,\n"
    code += "                       const PIKG::S32* __restrict__ ni,\n"
    code += "                       #{$force_name}**  __restrict__  force)\n"
    code += "{\n"
    code += "    prof.start(\"Retrieve#{$kernel_name}\");\n"
    code += "    prof.increaseLevel();\n"

    code += "    const int thread_id = 0;\n"
    code += "    for (size_t sw_off=0, i=0; sw_off< swalks.size(); sw_off+=BSIZE, i++){\n"
    code += "      FORCEBuffer& hFORCE = pdatalist[thread_id].hFORCE[i];\n"
    code += "      prof.start(\"ReadForce\");\n"
    code += "      status = clEnqueueReadBuffer(queue, dFORCE[i], CL_TRUE, 0, #{io_total_size[2]/32}*sizeof(PIKG::F32)*NIMAX*BSIZE, hFORCE.buffer, 1, &tail_kernel_event[i], nullptr); PIKG_MN_CHECK(status);\n"
    code += "      prof.end(\"ReadForce\");\n"
    code += "      prof.start(\"RetrieveBuffer\");\n"
    code += "      retrieve_transfer_buffer(i,swalks, sw_off, thread_id);\n"
    code += "      prof.end(\"RetrieveBuffer\");\n"
    code += "      prof.start(\"ReleaseMNCLBuffer\");\n"
    #fvars.each{ |v|
    #  iotype   = h[v][0]
    #  type     = h[v][1]
    #  fdpsname = h[v][2]
    #  if iotype =~ /EPJ/
    #    read_or_write = "CL_MEM_WRITE_ONLY"
    #    read_or_write = "CL_MEM_READ_ONLY" if iotype == "FORCE"
    #    size = "NIMAX*BSIZE"
    #    size = "NJMAX*BSIZE*SIZEOF_NJ_CLUSTER" if iotype == "EPJ"
    #    suffixes = [""]
    #    suffixes = ["B0","B1"] if iotype == "EPJ"
    #    get_vector_elements(type).each{ |s|
    #      suffixes.each{ |suffix|
    #        code += "      status = clReleaseMemObject(#{iotype}#{fdpsname}#{s}Buffer#{suffix}[i]); PIKG_MN_CHECK(status);\n"
    #      }
    #    }
    #  end
    #}
    code += "      status = clReleaseMemObject(dEPI[i]); PIKG_MN_CHECK(status);\n"
    code += "      status = clReleaseMemObject(dEPJB0[i]); PIKG_MN_CHECK(status);\n"
    code += "      status = clReleaseMemObject(dEPJB1[i]); PIKG_MN_CHECK(status);\n"
    code += "      status = clReleaseMemObject(dFORCE[i]); PIKG_MN_CHECK(status);\n"
    code += "      prof.end(\"ReleaseMNCLBuffer\");\n"
    code += "    }\n"

    code += "    int ni_max = 0;\n"
    code += "    for(int w=0;w<n_walk;w++) ni_max = std::max(ni_max,ni[w]);\n"
    code += "    #pragma omp parallel for\n"
    code += "    for (int index=0; index<n_walk*ni_max; index++){\n"
    code += "      int w = index / ni_max;\n"
    code += "      int i = index % ni_max;\n"
    code += "      if(i<ni[w]){\n"
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      if iotype == "FORCE"
        get_vector_elements(type).each{ |s|
          code += "             force[w][i].#{fdpsname}#{".#{s}" if s != ""} += store[w][i].#{fdpsname}#{".#{s}" if s != ""};\n"
        }
      end
    }
    code += "      }\n"
    code += "    }\n"

    code += "    prof.decreaseLevel();\n"
    code += "    prof.end(\"Retrieve#{$kernel_name}\");\n"
    code += "    return 0;\n"
    code += "}\n"

    code += "int make_transfer_buffer_index(const int tag, const std::vector<Subwalk> &swalks, const #{$epi_name}* epi[], const PIKG::S32* id_epj[], int sw_off, int thread_id ){\n"

    code += "  prof.increaseLevel();\n"

    code += "  //calculate size (NJ)\n"
    code += "  int nj_block_size_max = 0;\n"
    code += "  for (size_t sw=0; sw<BSIZE; sw++){\n"
    code += "    if (sw+sw_off >= swalks.size() ) break;\n"
    code += "    const Subwalk& walk = swalks[sw+sw_off];\n"
    code += "    nj_block_size_max = std::max(walk.nj_size, nj_block_size_max);\n"
    code += "  }\n"
    code += "  int SIZEOF_NJ_CLUSTER = (#{64/$min_element_size}*LW_PER_DAR_ENTRY*nj_block_size_max + NJMAX - 1) / NJMAX;\n"

    code += "  // reallocate memory if lack of memory size\n"
    code += "  if ( SIZEOF_NJ_CLUSTER > pdatalist[thread_id].NJ_CLUSTER_MAX[tag] ){\n"
    code += "   pdatalist[thread_id].reallocate_epj_index(tag,SIZEOF_NJ_CLUSTER * 2);\n"
    code += "  }\n"

    code += "  EPIBuffer&   hEPI   = pdatalist[thread_id].hEPI[tag];\n"
    code += "  PIKG::S32*   hIdEPJ = pdatalist[thread_id].hIdEPJ[tag];\n"
    code += "  FORCEBuffer& hFORCE = pdatalist[thread_id].hFORCE[tag];\n"

    code += "  prof.start(\"CopyEPI\");\n"
    code += "  #pragma omp parallel for\n"
    code += "  for (int i=0; i<NIMAX*BSIZE; i++){\n"
    code += "    const int sw = i/NIMAX;\n"
    code += "    if ( sw+sw_off < swalks.size() ){\n"
    code += "      const auto& walk = swalks[sw+sw_off];\n"
    code += "      const int w = walk.walk_id;\n"
    code += "      const int ni_off = walk.ni_off;\n"
    code += "      const int ii = i%NIMAX;\n"
    if $min_element_size == 32
      code += "      const int index32 = (BSIZE*(ii/(2*l2bm_chunk_size)) + sw)*(2*l2bm_chunk_size) + ((ii/2)%l2bm_chunk_size)*2 + (1^(ii%2));\n"
    end
    if $max_element_size == 64
      code += "      const int index64 = (BSIZE*(ii/l2bm_chunk_size) + sw)*l2bm_chunk_size + ii%l2bm_chunk_size;\n"
    end
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      size = get_single_data_size(type)
      if iotype == "EPI"
        get_vector_elements(type).each{ |s|
          code += "      h#{iotype}.#{fdpsname}#{s}[index#{size}] = epi[w][ni_off + ii].#{fdpsname}"
          code += ".#{s}" if s != ""
          code += ";\n"
        }
      end
    }
    code += "    }\n"
    code += "  }\n"
    code += "  prof.end(\"CopyEPI\");\n"

    code += "  prof.start(\"CopyIdEPJ\");\n"
    code += "  #pragma omp parallel for\n"
    code += "  for (int j=0; j<NJBLOCKMAX*BSIZE; j++){\n"
    code += "    const int sw = j/NJBLOCKMAX;\n"
    code += "    if(sw+sw_off < swalks.size()){\n"
    code += "      const auto& swalk = swalks[sw+sw_off];\n"
    code += "      const int w = swalk.walk_id;\n"
    code += "      const int jj = j%NJBLOCKMAX;\n"
    # EPJ is distributed by l1bmd. each pe has data in order of 64*j + $peid
    # DRAM indirect access is unit of 16 LW with 32 bit integer entry. 4 entries per L1B.
    code += "      int index = (4*BSIZE)*(jj/4) + 4*sw + 2*((jj/2)%2) + (1^(jj%2));\n"
    code += "      int offset = DRAM_J_INDIRECT_OFFSET/16;\n"
    code += "      const int increment = DRAM_INDIRECT_OFFSET_SIZE/16;\n"
    code += "      if(jj < swalk.nj_size){\n"
    for j in 1..nelem_epj do
      code += "        hIdEPJ[index] = offset + id_epj[w][jj];\n"
      code += "        index += NJBLOCKMAX*BSIZE;\n"
      code += "        offset += increment;\n"
    end
    code += "      }else{\n"
    for j in 1..nelem_epj do
      code += "        hIdEPJ[index] = offset + ID_EPJ_DUMMY;\n"
      code += "        index += NJBLOCKMAX*BSIZE;\n"
      code += "        offset += increment;\n"
    end
    code += "      }\n"
    code += "    }\n"
    code += "  }\n"
    code += "  prof.end(\"CopyIdEPJ\");\n"

    code += "  prof.decreaseLevel();\n"
    code += "  return SIZEOF_NJ_CLUSTER;\n"
    code += "}\n"


    code += "void DEVICE_RUN_DIVJ_INDEX(int tag, int Nloop, int NI, int NJpPE, int thread_id, const std::vector<Subwalk> &swalks, const PIKG::S32* id_epj[], int sw_off){ // Nloop (# of J accumeration)\n"
    code += "  EPIBuffer&   hEPI   = pdatalist[thread_id].hEPI[tag];\n"
    code += "  FORCEBuffer& hFORCE = pdatalist[thread_id].hFORCE[tag];\n"
    code += "  PIKG::S32*   hIdEPJ   = pdatalist[thread_id].hIdEPJ[tag];\n"

    code += "  prof.start(\"LoadKernels\");\n"
    code += "  cl_kernel Hkernel  = load_kernel(context,\"packs/nbY1_HI\"   + std::to_string(tag) + \".vsm\",\"Hkernel\"  + std::to_string(tag) + \"_index\");\n"
    code += "  cl_kernel FFkernel = load_kernel(context,\"packs/nbY1_MSI\"  + std::to_string(tag) + \".vsm\",\"FFkernel\" + std::to_string(tag) + \"_index\");\n"
    code += "  cl_kernel RRkernel = load_kernel(context,\"packs/nbY1_MSIR\" + std::to_string(tag) + \".vsm\",\"RRkernel\" + std::to_string(tag) + \"_index\");\n"
    code += "  cl_kernel FTkernel = load_kernel(context,\"packs/nbY1_MT\"  + std::to_string(tag) + \".vsm\",\"FTkernel\" + std::to_string(tag) + \"_index\");\n"
    code += "  cl_kernel RTkernel = load_kernel(context,\"packs/nbY1_MRT\" + std::to_string(tag) + \".vsm\",\"RTkernel\" + std::to_string(tag) + \"_index\");\n"
    code += "  prof.end(\"LoadKernels\");\n"

    code += "  // NI SET\n"
    code += "  prof.start(\"WriteEPI\");\n"
    code += "  status = clEnqueueWriteBuffer(queue, dEPI[tag], CL_FALSE, 0, #{io_total_size[0]/32}*sizeof(PIKG::F32)*NIMAX*BSIZE, hEPI.buffer, 0, nullptr, &write_events[tag][0]); PIKG_MN_CHECK(status);\n"
    code += "  prof.end(\"WriteEPI\");\n"

    code += "  // NJ SET\n"
    code += "  prof.start(\"WriteIdEPJFirst\");\n"
    code += "  status = clEnqueueWriteBuffer(queue, dEPJB0[tag], CL_FALSE, 0, #{nelem_epj}*sizeof(PIKG::F32)*NJBLOCKMAX*BSIZE, hIdEPJ, 0, nullptr, &write_events[tag][1]); PIKG_MN_CHECK(status);\n"
    code += "  prof.end(\"WriteIdEPJFirst\");\n"

    code += "  prof.start(\"EnqueueHeadKernel\");\n"
    code += "  status = clEnqueueTask(queue, Hkernel, 2, write_events[tag], &head_kernel_event[tag]); PIKG_MN_CHECK(status);\n"
    code += "  prof.end(\"EnqueueHeadKernel\");\n"
   
    code += "  assert(Nloop <= NLOOP_MAX);\n" 
    code += "  for ( int loop=0; loop< Nloop ; loop++ ){\n"

    code += "    // NJ SET\n"
    code += "    if (loop<Nloop-1){\n"
    code += "#ifdef PIKG_PROFILE_COPY_EP_TO_BUFFER\n"
    code += "      auto s = std::chrono::high_resolution_clock::now();\n"
    code += "#endif\n"
    code += "      prof.start(\"CopyIdEPJtoBuffer\");\n"
    code += "      #pragma omp parallel for\n"
    code += "      for (int j=0; j<NJBLOCKMAX*BSIZE; j++){\n"
    code += "        const int sw = j/NJBLOCKMAX;\n"
    code += "        if(sw+sw_off < swalks.size()){\n"
    code += "          const int jj = j%NJBLOCKMAX;\n"
    code += "          const auto swalk = swalks[sw+sw_off];\n"
    code += "          const int w = swalk.walk_id;\n"
    code += "          const int nj = swalk.nj_size;\n"
    code += "          const int jjj = (loop+1)*NJBLOCKMAX + jj;\n"
    code += "          int index = (loop+1)*#{nelem_epj}*NJBLOCKMAX*BSIZE + (4*BSIZE)*(jj/4) + 4*sw + 2*((jj/2)%2) + (1^(jj%2));\n"
    code += "          int offset = DRAM_J_INDIRECT_OFFSET/16;\n"
    code += "          const int increment = DRAM_INDIRECT_OFFSET_SIZE/16;\n"
    code += "          if(jjj < nj){\n"
    for j in 0..(nelem_epj-1) do
      code += "            hIdEPJ[index] = offset + id_epj[w][jjj];\n"
      code += "            index += NJBLOCKMAX*BSIZE;\n"
      code += "            offset += increment;\n"
    end
    code += "          }else{\n"
    for j in 0..(nelem_epj-1) do
      code += "            hIdEPJ[index] = offset + ID_EPJ_DUMMY;\n"
      code += "            index += NJBLOCKMAX*BSIZE;\n"
      code += "            offset += increment;\n"
    end
    code += "          }\n"
    code += "        }\n"
    code += "      }\n"

    code += "      prof.end(\"CopyIdEPJtoBuffer\");\n"
    code += "#ifdef PIKG_PROFILE_COPY_EP_TO_BUFFER\n"
    code += "      auto e = std::chrono::high_resolution_clock::now();\n"
    code += "      auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(e-s).count() * 1e-9; // in sec\n"
    code += "      std::cout << \"ElapsedTime(sec): \" << elapsed << \" Bandwidth(GB/s): \" << 2*NJMAX*BSIZE*#{epj_total_size}/elapsed * 1e-9 << std::endl;\n"
    code += "#endif\n"
    code += "      cl_event& wait_kernel_event = (loop==0) ? head_kernel_event[tag] : kernel_events[tag][loop-1];\n"
    code += "      prof.start(\"WriteEPJ\");\n"
    code += "      cl_mem& dIdEPJ = (loop%2==0) ? dEPJB1[tag] : dEPJB0[tag];\n"
    code += "      status = clEnqueueWriteBuffer(queue, dIdEPJ, CL_FALSE, 0, #{nelem_epj}*sizeof(PIKG::F32)*NJBLOCKMAX*BSIZE, hIdEPJ+#{nelem_epj}*NJBLOCKMAX*BSIZE*(loop+1), 1, &wait_kernel_event, &write_events[tag][1 + (loop+1)]); PIKG_MN_CHECK(status);\n"
    code += "      prof.end(\"WriteEPJ\");\n"
    code += "    }\n"
    code += "    prof.start(\"EnqueueMainKernel\");\n"
    code += "    if(loop==Nloop-1){\n"
    code += "      cl_kernel& kernel = (loop%2==0)? FTkernel : RTkernel;\n"
    code += "      if(loop == 0){\n"
    code += "        status = clEnqueueTask(queue, kernel, 1, &head_kernel_event[tag], &tail_kernel_event[tag]); PIKG_MN_CHECK(status);\n"
    code += "      }else{\n"
    code += "        status = clEnqueueTask(queue, kernel, 1, &kernel_events[tag][loop-1], &tail_kernel_event[tag]); PIKG_MN_CHECK(status);\n"
    code += "      }\n"
    code += "    }else{\n"
    code += "      cl_kernel& kernel = (loop%2==0)? FFkernel : RRkernel;\n"
    code += "      status = clEnqueueTask(queue, kernel, 1, &write_events[tag][1+(loop+1)], &kernel_events[tag][loop]); PIKG_MN_CHECK(status);\n"
    code += "    }\n"
    code += "    prof.end(\"EnqueueMainKernel\");\n"
    code += "  }\n"
    code += "  pdatalist[thread_id].reset_epj_adr(tag);\n"
    code += "}\n"


    code += "void RUN_WALKS_INDEX(\n"
    code += "                             #{$force_name}     *store[],\n"
    code += "                             const PIKG::S32      n_walk,\n"
    code += "                             const #{$epi_name} *epi[],\n"
    code += "                             const PIKG::S32      n_epi[],\n"
    code += "                             const PIKG::S32    *id_epj[],\n" # blocked id
    code += "                             const PIKG::S32      n_epj[]){\n" # number of blocks

    code += "    if ((NJMAX%#{64*(64/$min_element_size)})!=0) { std::cout << \"NJ != 0 mod #{64*(64/$min_element_size)} Error!\" << std::endl; exit(0); }\n"
    code += "    prof.start(\"RunWalks\");\n"
    code += "    prof.increaseLevel();\n"

    code += "    swalks.clear();\n"
    code += "    prof.start(\"DevideWalks\");\n"
    code += "    divide_walks( swalks, n_walk, n_epi, n_epj);\n"
    code += "    prof.end(\"DevideWalks\");\n"

    code += "    for (size_t sw_off=0, tag=0; sw_off< swalks.size(); sw_off+=BSIZE, tag++){\n"
    code += "      assert(tag < MAX_KERNEL_LAUNCH);\n"
    code += "      int thread_id = omp_get_thread_num();\n"
    code += "      assert(thread_id == 0);\n"

    code += "      prof.start(\"MakeTransferBuffer\");\n"
    code += "      int SIZEOF_NJ_CLUSTER = make_transfer_buffer_index(tag,swalks, epi, id_epj, sw_off, thread_id);\n"
    code += "      prof.end(\"MakeTransferBuffer\");\n"
    code += "      prof.start(\"CreateMNCLBuffer\");\n"
    code += "      {\n"
    code += "        dEPI[tag] = clCreateBufferWithAttributes(context, CL_MEM_READ_ONLY, 0, sizeof(PIKG::F64)*(tag*DRAM_KERNEL_OFFSET+DRAM_IN_OFFSET), #{io_total_size[0]/32}*sizeof(PIKG::F32)*NIMAX*BSIZE, nullptr, &status); PIKG_MN_CHECK(status);\n"
    code += "        dEPJB0[tag] = clCreateBufferWithAttributes(context, CL_MEM_READ_ONLY, 0, sizeof(PIKG::F64)*(tag*DRAM_KERNEL_OFFSET+DRAM_IN_J_OFFSET), #{nelem_epj}*sizeof(PIKG::F32)*NJMAX*BSIZE*SIZEOF_NJ_CLUSTER, nullptr, &status); PIKG_MN_CHECK(status);\n"
    code += "        dEPJB1[tag] = clCreateBufferWithAttributes(context, CL_MEM_READ_ONLY, 0, sizeof(PIKG::F64)*(tag*DRAM_KERNEL_OFFSET+DRAM_IN_J_OFFSET+DRAM_OFFSET_SIZE), #{nelem_epj}*sizeof(PIKG::F32)*NJMAX*BSIZE*SIZEOF_NJ_CLUSTER, nullptr, &status); PIKG_MN_CHECK(status);\n"
    code += "        dFORCE[tag] = clCreateBufferWithAttributes(context, CL_MEM_READ_ONLY, 0, sizeof(PIKG::F64)*(tag*DRAM_KERNEL_OFFSET+DRAM_OUT_OFFSET), #{io_total_size[2]/32}*sizeof(PIKG::F32)*NIMAX*BSIZE, nullptr, &status); PIKG_MN_CHECK(status);\n"
    code += "      }\n"
    code += "      prof.end(\"CreateMNCLBuffer\");\n"

    code += "      prof.start(\"DeviceRunDivJIndex\");\n"
    code += "      prof.increaseLevel();\n"
    code += "      DEVICE_RUN_DIVJ_INDEX( tag, SIZEOF_NJ_CLUSTER, NIMAX, NJpPE, thread_id, swalks, id_epj, sw_off);\n"
    code += "      prof.decreaseLevel();\n"
    code += "      prof.end(\"DeviceRunDivJIndex\");\n"

    code += "    }\n"
    code += "    prof.decreaseLevel();\n"
    code += "    prof.end(\"RunWalks\");\n"
    code += "}\n"

    code += "PIKG::S32 Dispatch#{$kernel_name}Index(\n"
    code += "                             const PIKG::S32          tag,\n"
    code += "                             const PIKG::S32          n_walk,\n"
    code += "                             const #{$epi_name}** __restrict__   epi,\n"
    code += "                             const PIKG::S32*    __restrict__   n_epi,\n"
    code += "                             const PIKG::S32** __restrict__   id_epj,\n"
    code += "                             const PIKG::S32*    __restrict__   n_epj,\n"
    code += "                             const #{$epj_name}* __restrict__   epj,\n"
    code += "                             const PIKG::S32   nsend_epj,\n"
    code += "                             const bool  send_flag){\n"
    code += "    prof.start(\"Dispatch#{$kernel_name}Index\");\n"
    code += "    prof.increaseLevel();\n"
    code += "    static bool is_initialized = false;\n"
    code += "    if(!is_initialized){\n"
    code += "      prof.start(\"Initialization\");\n"
    code += "      Initialize_this_header();\n"
    code += "      is_initialized = true;\n"
    code += "      prof.end(\"Initialization\");\n"
    code += "      // create epj buffer\n"
    code += "      int epj_offset   = sizeof(PIKG::F64)*(DRAM_J_INDIRECT_OFFSET);\n"
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      stype = get_single_element_type(type)
      if iotype =~ /EPJ/
        read_or_write = "CL_MEM_WRITE_ONLY"
        size = "sizeof(PIKG::F64)*DRAM_INDIRECT_OFFSET_SIZE"
        increment = "sizeof(PIKG::F64)*DRAM_INDIRECT_OFFSET_SIZE"
        suffixes = ["B0"]
        get_vector_elements(type).each{ |s|
          offset = "epj_offset"
          suffixes.each{ |suffix|
            isx = ""
            isx = " + sizeof(PIKG::F64)*DRAM_OFFSET_SIZE" if suffix == "B1"
            code += "      #{iotype}#{fdpsname}#{s}Buffer#{suffix}[tag] = clCreateBufferWithAttributes(context, #{read_or_write}, 0, #{offset}#{isx}, #{size}, nullptr, &status); PIKG_MN_CHECK(status);\n"
          }
          code += "      #{offset} += #{increment};\n"
        }
      end
    }
    code += "    }\n"
    code += "    if(send_flag){\n"
    code += "      prof.start(\"SendEPJ\");\n"
    code += "      const int nsend_max = ((nsend_epj+#{2*16*(64/$min_element_size)-1})/#{16*(64/$min_element_size)})*#{16*(64/$min_element_size)};\n" # dummy particles must be added
    code += "      assert(nsend_max/#{64/$max_element_size} < DRAM_INDIRECT_OFFSET_SIZE);\n"
    code += "      pdatalist[0].reallocate_epj_for_indirect(nsend_max);\n"
    code += "      int epj_offset = sizeof(PIKG::F64)*DRAM_J_INDIRECT_OFFSET;\n"

    code += "      // copy epj to local buffer\n"
    code += "      EPJBuffer&   hEPJ   = pdatalist[0].hEPJ[0];\n"
    code += "      const #{$epj_name} epj_dummy;\n"
    code += "      #pragma omp parallel for\n"
    code += "      for(int j=0;j<nsend_max;j++){\n"
    code += "        if(j < nsend_epj){\n"
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      if iotype == "EPJ"
        get_vector_elements(type).each{ |s|
          code += "            h#{iotype}.#{fdpsname}#{s}[j] = epj[j].#{fdpsname}"
          code += ".#{s}" if s != ""
          code += ";\n"
        }
      end
    }
    code += "        }else{\n"
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      if iotype == "EPJ"
        get_vector_elements(type).each{ |s|
          code += "            h#{iotype}.#{fdpsname}#{s}[j] = epj_dummy.#{fdpsname}"
          code += ".#{s}" if s != ""
          code += ";\n"
        }
      end
    }
    code += "        }\n"
    code += "      }\n"
    write_count = 0
    fvars.each{ |v|
      iotype   = h[v][0]
      type     = h[v][1]
      fdpsname = h[v][2]
      stype = get_single_element_type(type)
      if iotype =~ /EPJ/
        read_or_write = "CL_MEM_WRITE_ONLY"
        increment = "sizeof(PIKG::F64)*DRAM_J_INDIRECT_OFFSET"
        suffixes = ["B0"]
        get_vector_elements(type).each{ |s|
          offset = "epj_offset"
          suffixes.each{ |suffix|
            code += "      status = clEnqueueWriteBuffer(queue, #{iotype}#{fdpsname}#{s}Buffer#{suffix}[tag], CL_FALSE, 0, sizeof(PIKG::#{stype})*nsend_max, hEPJ.#{fdpsname}#{s}, 0, nullptr, &write_events[tag][#{write_count}]); PIKG_MN_CHECK(status);\n"
            write_count+=1
          }
          code += "      #{offset} += #{increment};\n"
        }
      end
    }
    code += "      waitForEventsByDummyRead(queue,#{write_count},write_events[tag]);\n"
    code += "      ID_EPJ_DUMMY = (nsend_max-#{16*(64/$min_element_size)})/#{16*(64/$min_element_size)};\n"
    code += "      prof.end(\"SendEPJ\");\n"
    code += "      prof.decreaseLevel();\n"
    code += "      prof.end(\"Dispatch#{$kernel_name}Index\");\n"
    code += "      return 0;\n"
    code += "    }\n"
    code += "    prof.start(\"PrepareStoreBuffer\");\n"
    code += "    assert(n_walk<=N_WALK_LIMIT);\n"
    code += "    assert(n_epi[0]>0 && n_epj[0]>0);\n"
    code += "    static bool first = true;\n"
    code += "    if(first){\n"
    code += "      store = new #{$force_name}*[N_WALK_LIMIT];\n"
    code += "      for (int iw=0; iw<N_WALK_LIMIT; iw++){\n"
    code += "        store[iw] = new #{$force_name}[N_EPI_MAX];\n"
    code += "      }\n"
    code += "      first = false;\n"
    code += "    }\n"
    code += "    prof.end(\"PrepareStoreBuffer\");\n"

    code += "    prof.start(\"CompressIndex\");\n"
    code += "    // sort and compress id_epj\n"
    code += "    PIKG::S32** id_epj_comp = new PIKG::S32*[n_walk];\n"
    code += "    PIKG::S32 n_epj_max = 0;\n"
    code += "    for (int w=0; w<n_walk; w++){\n"
    code += "      id_epj_comp[w] = new PIKG::S32[n_epj[w]];\n"
    code += "      n_epj_max = std::max(n_epj_max,n_epj[w]);\n"
    code += "    }\n"
    code += "    #pragma omp parallel for\n"
    code += "    for (int j=0; j<n_epj_max*n_walk; j++){\n"
    code += "      const int w = j/n_epj_max;\n"
    code += "      const int jj = j%n_epj_max;\n"
    code += "      const int nj = n_epj[w];\n"
    code += "      if(jj < nj){\n"
    code += "        id_epj_comp[w][jj] = id_epj[w][jj] / #{16*(64/$min_element_size)};\n"
    code += "      }\n"
    code += "    }\n"
    code += "    PIKG::S32 *n_epj_block;\n n_epj_block = new PIKG::S32[n_walk];\n"
    code += "    #pragma omp parallel for\n"
    code += "    for(int w=0;w<n_walk;w++){\n"
    code += "      std::sort(id_epj_comp[w],id_epj_comp[w]+n_epj[w]);\n"
    code += "      auto last = std::unique(id_epj_comp[w],id_epj_comp[w]+n_epj[w]);\n"
    code += "      n_epj_block[w] = (PIKG::S32)(last - id_epj_comp[w]);\n"
    code += "    }\n"
    code += "    prof.end(\"CompressIndex\");\n"

    code += "    RUN_WALKS_INDEX( store, n_walk, epi, n_epi, (const PIKG::S32**)id_epj_comp, n_epj_block);\n"
    code += "    prof.decreaseLevel();\n"
    code += "    prof.end(\"Dispatch#{$kernel_name}Index\");\n"

    code += "    for(int w=0;w<n_walk;w++){\n"
    code += "      delete[] id_epj_comp[w];\n"
    code += "    }\n"
    code += "    delete[] id_epj_comp;\n"
    code += "    delete[] n_epj_block;\n"
    code += "    return 0;\n"
    code += "}\n"

    code += "struct #{$kernel_name}{\n"
    $varhash.each_with_index{|v,i|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        code += "  static " + get_declare_type(type,conversion_type) + " " + name + ";\n"
      end
    }
    code += "  #{$kernel_name}("
    needComma = nil
    $varhash.each_with_index{|v,i|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        code += "," if needComma
        needComma = 1
        code += "const " + get_declare_type(type,conversion_type) + " _" + name
      end
    }
    code += "){\n"
    code += "    initialize_#{$kernel_name}("
    needComma = nil
    $varhash.each_with_index{|v,i|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        code += "," if needComma
        needComma = 1
        code += "_" + name
      end
    }
    code += ");\n"
    code += "  }\n"
    code += "  ~#{$kernel_name}(){\n"
    code += "#ifdef PIKG_ENABLE_PERFETTO_TRACE\n"
    code += "    grape_pfn::base::StopTrace();\n"
    code += "#endif\n"
    code += "  }\n"
    code += "  void initialize_#{$kernel_name}("
    count = 0
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        code += "," if count > 0
        code += "const " + get_declare_type(type,conversion_type) + " " + name + "_"
        count += 1
      end
    }
    code += "){\n"
    code += "      // do nothing here because dynamic setting of member variable is not allowed with PIKG now\n"
    code += "  }\n"
    code += "  void operator()(const #{$epi_name}* epi,\n"
    code += "                  const int nepi,\n"
    code += "                  const #{$epj_name}* epj,\n"
    code += "                  const int nepj,\n"
    code += "                  #{$force_name}* force){\n"
    code += "    prof.start(\"Total\");\n"
    code += "    prof.increaseLevel();\n"

    code += "    Dispatch#{$kernel_name}(0,1,&epi,&nepi,&epj,&nepj);\n"
    code += "    Retrieve#{$kernel_name}(0,1,&nepi,&force);\n"

    code += "    prof.decreaseLevel();\n"
    code += "    prof.end(\"Total\");\n"
    code += "  }\n"
    code += "};\n"

    code += "void clearTimeProfiler(){ prof.clearAll(); }\n"

    File.open($output_file, mode = 'w'){ |f|
      f.write(code)
    }
  end
end

