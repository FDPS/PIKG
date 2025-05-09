class Kernelprogram
  def reserved_func_def_cuda(conversion_type)
    abort if conversion_type != "CUDA"
    code = ""
    code += "__host__ __device__ PIKG::F64 inv(PIKG::F64 op){ return 1.0/op; }\n"
    code += "__host__ __device__ PIKG::F32 inv(PIKG::F32 op){ return 1.f/op; }\n"

    code += "__host__ __device__ PIKG::F64 table(PIKG::F64 tab[],PIKG::S64 i){ return tab[i]; }\n"
    code += "__host__ __device__ PIKG::F32 table(PIKG::F32 tab[],PIKG::S32 i){ return tab[i]; }\n"

    code += "__host__ __device__ PIKG::F64 to_float(PIKG::U64 op){return (PIKG::F64)op;}\n"
    code += "__host__ __device__ PIKG::F32 to_float(PIKG::U32 op){return (PIKG::F32)op;}\n"
    code += "__host__ __device__ PIKG::F64 to_float(PIKG::S64 op){return (PIKG::F64)op;}\n"
    code += "__host__ __device__ PIKG::F32 to_float(PIKG::S32 op){return (PIKG::F32)op;}\n"
    code += "__host__ __device__ PIKG::S64   to_int(PIKG::F64 op){return (PIKG::S64)op;}\n"
    code += "__host__ __device__ PIKG::S32   to_int(PIKG::F32 op){return (PIKG::S32)op;}\n"
    code += "__host__ __device__ PIKG::U64  to_uint(PIKG::F64 op){return (PIKG::U64)op;}\n"
    code += "__host__ __device__ PIKG::U32  to_uint(PIKG::F32 op){return (PIKG::U32)op;}\n"

    code += "template<typename T> __host__ __device__ PIKG::F64 to_f64(const T& op){return (PIKG::F64)op;}\n"
    code += "template<typename T> __host__ __device__ PIKG::F32 to_f32(const T& op){return (PIKG::F32)op;}\n"
    #code += "template<typename T> __host__ __device__ PIKG::F16 to_f16(const T& op){return (PIKG::F16)op;}\n"
    code += "template<typename T> __host__ __device__ PIKG::S64 to_s64(const T& op){return (PIKG::S64)op;}\n"
    code += "template<typename T> __host__ __device__ PIKG::S32 to_s32(const T& op){return (PIKG::S32)op;}\n"
    #code += "template<typename T> __host__ __device__ PIKG::S16 to_s16(const T& op){return (PIKG::S16)op;}\n"
    code += "template<typename T> __host__ __device__ PIKG::U64 to_u64(const T& op){return (PIKG::U64)op;}\n"
    code += "template<typename T> __host__ __device__ PIKG::U32 to_u32(const T& op){return (PIKG::U32)op;}\n"
    #code += "template<typename T> __host__ __device__ PIKG::U16 to_u16(const T& op){return (PIKG::U16)op;}\n"

    code
  end

  def generate_optimized_cuda_kernel(conversion_type,h = $varhash)
    abort "error: --class-file option must be specified for conversion_type CUDA" if $epi_file == nil
    
    nepjsh = 1000
    $max_element_size = 1
    $min_element_size = 1
    
    code = String.new

    code += "#include <cuda/std/limits>\n"
    code += "#include \"#{$epi_file}\"\n"
    code += "#include \"#{$epj_file}\"\n" if $epj_file != $epi_file
    code += "#include \"#{$force_file}\"\n" if ($force_file != $epi_file) && ($force_file != $epj_file)
    code += "#include \"pikg_cuda_pointer.hpp\"\n"
    code += "#include \"pikg_vector.hpp\"\n"
    code += "#include \"pikg_time_profiler.hpp\"\n"

    code += "static PIKG::TimeProfiler prof;\n"

    fvars = generate_force_related_map(@statements)
    # GPU class definition
    code += "struct EpiGPU{\n"
    fvars.each{ |v|
      iotype = h[v][0]
      modifier = h[v][3]
      if iotype == "EPI"
        type = h[v][1]
        modifier = h[v][3]
        if modifier == "local"
          fdpsname = v
        else
          fdpsname = h[v][2]
        end
        code += Declaration.new([type,fdpsname]).convert_to_code(conversion_type)
      end
    }
    code += "};\n"
    code += "struct EpjGPU{\n"
    fvars.each{ |v|
      iotype = h[v][0]
      modifier = h[v][3]
      if iotype == "EPJ"
        type = h[v][1]
        modifier = h[v][3]
        if modifier == "local"
          fdpsname = v
        else
          fdpsname = h[v][2]
        end
        code += Declaration.new([type,fdpsname]).convert_to_code(conversion_type)
      end
    }
    code += "};\n"
    code += "struct ForceGPU{\n"
    fvars.each{ |v|
      iotype = h[v][0]
      if iotype == "FORCE"
        type = h[v][1]
        modifier = h[v][3]
        if modifier == "local"
          fdpsname = v
        else
          fdpsname = h[v][2]
        end
        code += Declaration.new([type,fdpsname]).convert_to_code(conversion_type)
      end
    }
    code += "};\n"


    # constant definition
    code += "enum{\n"
    code += "  N_THREAD_GPU = 32,\n"
    code += "  N_WALK_LIMIT = 2048,\n"
    code += "  NI_LIMIT     = N_WALK_LIMIT*1000,\n"
    code += "  NJ_LIMIT     = N_WALK_LIMIT*100000,\n"
    code += "};\n"

    # h2d copy func
    #code += "void CopyEpToDevice(const Epi* epi,\n"
    #code += "                    const PIKG::S32 ni,\n"
    #code += "                    const #{$epj_name}* epj,\n"
    #code += "                    const PIKG::S32 nj){\n"
    #code += "  epi_dev.resize(ni);\n"
    #code += "  epi_dev.h2d(epi);\n"
    #code += "  epj_dev.resize(nj);\n"
    #code += "  epj_dev.h2d(epj);\n"
    #code += "  force_dev.resize(ni);\n"
    #code += "}\n"
    # d2h copy func
    #code += "void CopyForceFromDevice(const #{$force_name}* force,\n"
    #code += "                         const PIKG::S32 ni){\n"
    #code += "  force_dev.h2d(force);\n"
    #code += "}\n"

    member_vars = String.new
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        member_vars += "," + name
      end
    }
    member_decls = String.new
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        member_decls += ",\n"
        member_decls += get_declare_type(type,conversion_type) + " " + name
      end
    }
    code += reserved_func_def_cuda(conversion_type)

    split_vars = Array.new
    ["EPI","EPJ","FORCE"].each{ |io|
      fvars.each{ |v|
        iotype = h[v][0]
        if iotype == io
          type = h[v][1]
          fdpsname = h[v][2]
          modifier = h[v][3]
          n = get_single_data_size(type) / $min_element_size
          split_vars += [v] if n > 1
        end
      }
    }

    code += "inline __device__ ForceGPU inner_kernel(\n"
    code += "				     EpiGPU epi,\n"
    code += "				     EpjGPU epj,\n"
    code += "				     ForceGPU force"
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        code += ",\n				     "
        code += get_declare_type(type,conversion_type) + " " + name
      end
    }
    code += ")\n"
    code += "{\n"
    # describe most inner loop
    new_s = Array.new

    @statements.each{ |s|
      code += s.convert_to_code(conversion_type) + "\n" if s.class == TableDecl
      next if s.class == Pragma or s.class == TableDecl
      next if s.class == Statement and h[get_name(s)][3] == "local"
      new_s.push(s)
    }

    generate_jloop_body(new_s,fvars,split_vars,conversion_type).each{ |s|
      tmp_s = s
      if get_name(tmp_s) != nil
        fvars.each{ |v|
          iotype = h[v][0]
          index = ["EPI","EPJ","FORCE"].index(iotype)
          if index
            type = h[v][1]
            fdpsname = h[v][2]
            modifier = h[v][3]
            if modifier == "local"
              replaced = ["epi","epj","force"][index] + "." + v
            else
              replaced = ["epi","epj","force"][index] + "." + fdpsname
            end
            name = tmp_s.name.replace_recursive(v,replaced)
            if s.class == Statement
              exp  = tmp_s.expression.replace_recursive(v,replaced)
              tmp_s = Statement.new([name,exp,s.type,s.op])
            end
          end
        }
      end
      code += tmp_s.convert_to_code(conversion_type) + "\n"
    }
    code += "  return force;\n"
    code += "}\n"
  
    code += "__device__ ForceGPU ForceKernel_1walk(\n"
    code += "				    EpjGPU *jpsh,\n"
    code += "				    const EpiGPU ipos,\n"
    code += "				    const int     id_walk,\n"
    code += "				    const int    *ij_disp,\n"
    code += "				    const EpjGPU *epj, \n"
    code += "				    ForceGPU accp"
    code += member_decls
    code += "){\n"
    code += "  const int tid = threadIdx.x;\n"
    code += "  const int j_head = ij_disp[id_walk  ];\n"
    code += "  const int j_tail = ij_disp[id_walk+1];\n"
    code += "  for(int j=j_head; j<j_tail; j+=N_THREAD_GPU){\n"
    code += "    jpsh[tid] = epj[j+tid];\n"
    code += "    if(j_tail-j < N_THREAD_GPU){\n"
    code += "      for(int jj=0; jj<j_tail-j; jj++){\n"
    code += "        accp = inner_kernel( ipos, jpsh[jj], accp#{member_vars});\n"
    code += "      }\n"
    code += "    }else{\n"
    code += "#pragma unroll 32\n"
    code += "      for(int jj=0; jj<N_THREAD_GPU; jj++){\n"
    code += "        accp = inner_kernel( ipos, jpsh[jj], accp#{member_vars});\n"
    code += "      }\n"
    code += "    }\n"
    code += "  }\n"
    code += "  return accp;\n"
    code += "}\n"

    code += "__device__ ForceGPU ForceKernel_2walk(\n"
    code += "				    EpjGPU jpsh[2][N_THREAD_GPU],\n"
    code += "				    const EpiGPU  ipos,\n"
    code += "				    const int     id_walk,\n"
    code += "				    const int     iwalk0,\n"
    code += "				    const int     iwalk1,\n"
    code += "				    const int    *ij_disp,\n"
    code += "				    const EpjGPU *epj, \n"
    code += "				    ForceGPU accp"
    code += member_decls
    code += "){\n"
    code += "  const int jbeg0 = ij_disp[iwalk0];\n"
    code += "  const int jbeg1 = ij_disp[iwalk1];\n"
    code += "  const int jend0 = ij_disp[iwalk0 + 1];\n"
    code += "  const int jend1 = ij_disp[iwalk1 + 1];\n"
    code += "  const int nj0   = jend0 - jbeg0;\n"
    code += "  const int nj1   = jend1 - jbeg1;\n"
    code += "  const int nj_longer  = nj0 > nj1 ? nj0 : nj1;\n"
    code += "  const int nj_shorter = nj0 > nj1 ? nj1 : nj0;\n"
    code += "  const int walk_longer= nj0 > nj1 ? 0 : 1;\n"
    code += "  const int jbeg_longer = nj0 > nj1 ? jbeg0 : jbeg1;\n"
    code += "  const int mywalk = id_walk==iwalk0 ? 0 : 1;\n"
    code += "  const int tid = threadIdx.x;\n"
    code += "  for(int j=0; j<nj_shorter; j+=N_THREAD_GPU){\n"
    code += "    jpsh[0][tid] = epj[jbeg0 + j + tid];\n"
    code += "    jpsh[1][tid] = epj[jbeg1 + j + tid];\n"
    code += "    if(nj_shorter-j < N_THREAD_GPU){\n"
    code += "      for(int jj=0; jj<nj_shorter-j; jj++){\n"
    code += "        accp = inner_kernel( ipos, jpsh[mywalk][jj], accp#{member_vars});\n"
    code += "      }\n"
    code += "    }else{\n"
    code += "#pragma unroll 32\n"
    code += "      for(int jj=0; jj<N_THREAD_GPU; jj++){\n"
    code += "        accp = inner_kernel( ipos, jpsh[mywalk][jj], accp#{member_vars});\n"
    code += "      }\n"
    code += "    }\n"
    code += "  }\n"
    code += "  for(int j=nj_shorter; j<nj_longer; j+=N_THREAD_GPU){\n"
    code += "    jpsh[0][tid] = epj[jbeg_longer + j + tid];\n"
    code += "    int jrem = nj_longer - j;\n"
    code += "    if(jrem < N_THREAD_GPU){\n"
    code += "      for(int jj=0; jj<jrem; jj++){\n"
    code += "	     if(mywalk == walk_longer)\n"
    code += "          accp = inner_kernel( ipos, jpsh[0][jj], accp#{member_vars});\n"
    code += "      }\n"
    code += "    }else{\n"
    code += "#pragma unroll 32\n"
    code += "      for(int jj=0; jj<N_THREAD_GPU; jj++){\n"
    code += "	     if(mywalk == walk_longer)\n"
    code += "          accp = inner_kernel( ipos, jpsh[0][jj], accp#{member_vars});\n"
    code += "      }\n"
    code += "    }\n"
    code += "  }\n"
    code += "  return accp;\n"
    code += "}\n"

    code += "__device__ ForceGPU ForceKernel_multiwalk(\n"
    code += "					const EpiGPU ipos,\n"
    code += "					const int     id_walk,\n"
    code += "					const int    *ij_disp,\n"
    code += "					const EpjGPU *epj, \n"
    code += "					ForceGPU accp"
    code += member_decls
    code += "){\n"
    code += "  const int j_head = ij_disp[id_walk  ];\n"
    code += "  const int j_tail = ij_disp[id_walk+1];\n"
    code += "\n"
    code += "  for(int j=j_head; j<j_tail; j++){\n"
    code += "    EpjGPU jp = epj[j];\n"
    code += "    accp = inner_kernel( ipos, jp, accp#{member_vars});\n"
    code += "  }\n"
    code += "  return accp;\n"
    code += "}\n"

    code += "__global__ void #{$kernel_name}_cuda(\n"
    code += "                  const int    * ij_disp,\n"
    code += "                  const int    * walk,\n"
    code += "                  const EpiGPU * epi,\n"
    code += "                  const EpjGPU * epj, \n"
    code += "                  ForceGPU     * force"
    code += member_decls
    code += ", bool clear = true"
    code +="){\n"
    code += "  int tid = blockDim.x * blockIdx.x + threadIdx.x;\n"
    # load epi
    code += "  const EpiGPU ip = epi[tid];\n"
    code += "  const int id_walk = walk[tid];\n"
    # init force
    code += "  ForceGPU accp;\n"
    fvars.each{ |v|
      iotype = h[v][0]
      if iotype == "FORCE"
        type = h[v][1]
        fdpsname = h[v][2]
        type_single = get_single_element_type(type)
        get_vector_elements(type).each_with_index{ |dim,j|
          dest = "accp." + fdpsname
          dest = Expression.new([:dot,dest,dim,type_single]) if dim != ""
          op = $accumhash[v][j]
          abort "accum_hash == nil" if $accumhash[v][j] == nil
          code += "    " + Duplicate.new([dest,get_initial_value(op,type_single,conversion_type),type_single]).convert_to_code("reference") + "\n"
        }
      end
    }
    code += "\n"
    code += "  if(clear){\n"
    fvars.each{ |v|
      iotype = h[v][0]
      if iotype == "FORCE"
        type = h[v][1]
        fdpsname = h[v][2]
        type_single = get_single_element_type(type)
        get_vector_elements(type).each_with_index{ |dim,j|
          dest = "force[tid]." + fdpsname
          dest = dest + "." + dim if dim != ""
          op = $accumhash[v][j]
          abort "accum_hash == nil" if $accumhash[v][j] == nil
          code += "    " + dest + "= 0.f;\n"
        }
      end
    }
    code += "  }\n"
    code += "  int t_head = blockDim.x * blockIdx.x;\n"
    code += "  int t_tail = t_head + N_THREAD_GPU - 1;\n"
    code += "  int nwalk_in_block = 1 + (walk[t_tail] - walk[t_head]);\n"
    code += "\n"
    code += "  __shared__ EpjGPU jpsh[2][N_THREAD_GPU];\n"
    code += "\n"
    code += "  if(1 == nwalk_in_block){\n"
    code += "    accp = ForceKernel_1walk(jpsh[0], ip, id_walk, ij_disp, epj, accp#{member_vars});\n"
    code += "  } else if(2 == nwalk_in_block){\n"
    code += "    int iwalk0 = walk[t_head];\n"
    code += "    int iwalk1 = walk[t_tail];\n"
    code += "    accp = ForceKernel_2walk(jpsh, ip, id_walk, iwalk0, iwalk1, ij_disp, epj, accp#{member_vars});\n"
    code += "  } else{\n"
    code += "    accp = ForceKernel_multiwalk(ip, id_walk, ij_disp, epj, accp#{member_vars});\n"
    code += "  }\n"
    # accumulate force
    fvars.each{ |v|
      iotype = h[v][0]
      if iotype == "FORCE"
        type = h[v][1]
        fdpsname = h[v][2]
        type_single = get_single_element_type(type)
        get_vector_elements(type).each_with_index{ |dim,j|
          dest = "force[tid]." + fdpsname
          suffix = ""
          suffix = "." + dim if dim != ""
          op = $accumhash[v][j]
          abort "accum_hash == nil" if $accumhash[v][j] == nil
          code += "    #{dest}#{suffix} += accp.#{fdpsname}#{suffix};\n"
        }
      end
    }
    code += "}\n"

    code += "static PIKG_CUDA::CUDAPointer<EpiGPU>   dev_epi;\n"
    code += "static PIKG_CUDA::CUDAPointer<EpjGPU>   dev_epj;\n"
    code += "static PIKG_CUDA::CUDAPointer<ForceGPU> dev_force;\n"
    code += "static PIKG_CUDA::CUDAPointer<int>  ij_disp;\n"
    code += "static PIKG_CUDA::CUDAPointer<int>   walk;\n"



    $varhash.each_with_index{|v,i|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        code += "static " + get_declare_type(type,conversion_type) + " " + name + ";\n"
      end
    }
    code += "void initialize_#{$kernel_name}("
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
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        type = v[1][1]
        code += name + "=" + name + "_;\n"
      end
    }
    code += "}\n"

    
    code += "PIKG::S32 Dispatch#{$kernel_name}(const PIKG::S32          tag,\n"
    code += "                   const PIKG::S32          n_walk,\n"
    code += "                   const #{$epi_name}    *epi[],\n"
    code += "                   const PIKG::S32          n_epi[],\n"
    code += "                   const #{$epj_name}    *epj[],\n"
    code += "                   const PIKG::S32          n_epj[],\n"
    if $spj_name != nil
      code += "               const #{$spj_name} *spj[],\n"
      code += "               const PIKG::S32        n_spj[],\n"
    end
    code += "                   bool clear = true)\n"
    code += "{\n"
    code += "  assert(n_walk <= N_WALK_LIMIT);\n"
    code += "  static bool init_call = true;\n"
    code += "  if(init_call){\n"
    code += "    dev_epi  .allocate(NI_LIMIT);\n"
    code += "    dev_epj  .allocate(NJ_LIMIT);\n"
    code += "    dev_force.allocate(NI_LIMIT);\n"
    code += "    ij_disp  .allocate(N_WALK_LIMIT+2);\n"
    code += "    walk     .allocate(NI_LIMIT);\n"
    code += "    init_call = false;\n"
    code += "  }\n"
    code += "  ij_disp[0] = 0;\n"
    code += "  for(int k=0; k<n_walk; k++){\n"
    code += "    ij_disp[k+1]  = ij_disp[k] + n_epj[k];\n"
    code += "    ij_disp[k+1] += n_spj[k];\n" if $spj_name != nil
    code += "  }\n"
    code += "  ij_disp[n_walk+1] = ij_disp[n_walk];\n"
    code += "\n"
    code += "  assert(ij_disp[n_walk] < NJ_LIMIT);\n"
    code += "  ij_disp.h2d(n_walk + 2);\n"
    code += "\n"
    code += "\n"
    code += "  int ni_tot = 0;\n"
    code += "  int nj_tot = 0;\n"
    code += "  prof.start(\"CopyEP\");\n"
    code += "  for(int iw=0; iw<n_walk; iw++){\n"
    code += "    for(int i=0; i<n_epi[iw]; i++){\n"
    # epi copy to epi_gpu
    fvars.each{ |v|
      iotype = h[v][0]
      modifier = h[v][3]
      if iotype == "EPI"
        type = h[v][1]
        modifier = h[v][3]
        if modifier == "local"
          fdpsname = v
          @statements.each{ |s|
            if get_name(s) == v
              dim = get_tail(s)
              new_name = "dev_epi[ni_tot]."+fdpsname
              new_name = new_name + "." + dim if dim != nil
              new_exp = s.expression.replace_fdpsname_recursive(h,true)
              code += "      " + Statement.new([new_name, new_exp,type]).convert_to_code("reference") + "\n"
            end
          }
        else
          fdpsname = h[v][2]
          get_vector_elements(type).each{|dim|
            if dim == ""
              code += "      " + Statement.new(["dev_epi[ni_tot]."+fdpsname,"epi[iw][i]."+fdpsname,type]).convert_to_code("reference") + "\n"
            else
              code += "      " + Statement.new(["dev_epi[ni_tot]."+fdpsname+"."+dim,"epi[iw][i]."+fdpsname+"."+dim,type]).convert_to_code("reference") + "\n"
            end
          }
        end
      end
    }
    code += "      walk[ni_tot] = iw;\n"
    code += "      ni_tot++;\n"
    code += "    }\n"
    code += "    for(int j=0; j<n_epj[iw]; j++){\n"
    # epj copy to epj_gpu
    fvars.each{ |v|
      iotype = h[v][0]
      modifier = h[v][3]
      if iotype == "EPJ"
        type = h[v][1]
        modifier = h[v][3]
        if modifier == "local"
          fdpsname = v
          @statements.each{ |s|
            if get_name(s) == v
              dim = get_tail(s)
              new_name = "dev_epj[nj_tot]."+fdpsname
              new_name = new_name + "." + dim if dim != nil
              new_exp = s.expression.replace_fdpsname_recursive(h,true)
              code += "      " + Statement.new([new_name, new_exp,type]).convert_to_code("reference") + "\n"
            end
          }
        else
          fdpsname = h[v][2]
          get_vector_elements(type).each{|dim|
            if dim == ""
              code += "      " + Statement.new(["dev_epj[nj_tot]."+fdpsname,"epj[iw][j]."+fdpsname,type]).convert_to_code("reference") + "\n"
            else
              code += "      " + Statement.new(["dev_epj[nj_tot]."+fdpsname+"."+dim,"epj[iw][j]."+fdpsname+"."+dim,type]).convert_to_code("reference") + "\n"
            end
          }
        end
      end
    }
    code += "      nj_tot++;\n"
    code += "    }\n"
    if $spj_name != nil
      code += "    for(int j=0; j<n_spj[iw]; j++){\n"
      # spj copy to spj_gpu
      fvars.each{ |v|
        iotype = h[v][0]
        if iotype == "EPJ"
          type = h[v][1]
          modifier = h[v][3]
          if modifier == "local"
            fdpsname = v
          else
            fdpsname = h[v][2]
          end
          get_vector_elements(type).each{|dim|
            if dim == ""
              code += "      " + Statement.new(["dev_epj[nj_tot]."+fdpsname,"spj[iw][j]."+fdpsname,type]).convert_to_code("reference") + "\n"
            else
              code += "      " + Statement.new(["dev_epj[nj_tot]."+fdpsname+"."+dim,"spj[iw][j]."+fdpsname+"."+dim,type]).convert_to_code("reference") + "\n"
            end
          }
        end
      }
      code += "      nj_tot++;\n"
      code += "    }\n"
    end
    code += "  }\n"
    code += "  prof.end(\"CopyEP\");\n"
    code += "  assert(ni_tot < NI_LIMIT);\n"
    code += "  int ni_tot_reg = ni_tot;\n"
    code += "  if(ni_tot_reg % N_THREAD_GPU){\n"
    code += "    ni_tot_reg /= N_THREAD_GPU;\n"
    code += "    ni_tot_reg++;\n"
    code += "    ni_tot_reg *= N_THREAD_GPU;\n"
    code += "  }\n"
    code += "  for(int i=ni_tot; i<ni_tot_reg; i++){\n"
    code += "    walk[i] = n_walk;\n"
    code += "  }\n"
    code += "\n"
    code += "  prof.start(\"MemcpyH2D\");\n"
    code += "  walk.h2d(ni_tot_reg);\n"
    code += "  dev_epi.h2d(ni_tot_reg);\n"
    code += "  dev_epj.h2d(nj_tot);\n"
    code += "  prof.end(\"MemcpyH2D\");\n"
    code += "\n"
    code += "  int nblocks  = ni_tot_reg / N_THREAD_GPU;\n"
    code += "  int nthreads = N_THREAD_GPU;\n"
    code += "#ifdef PIKG_MEASURE_CUDA_KERNEL_TIME\n"
    code += "  prof.start(\"Kernel\");\n"
    code += "#endif\n"
    code += "  #{$kernel_name}_cuda <<<nblocks, nthreads>>> (ij_disp, walk,  dev_epi, dev_epj, dev_force"
    $varhash.each{|v|
      iotype = v[1][0]
      if iotype == "MEMBER"
        name = v[0]
        code += ",#{name}"
      end
    }
    code += ",clear);\n"
    code += "#ifdef PIKG_MEASURE_CUDA_KERNEL_TIME\n"
    code += "  cudaDeviceSynchronize();\n"
    code += "  prof.end(\"Kernel\");\n"
    code += "#endif\n"
    code += "\n"
    code += "  return 0;\n"
    code += "}\n"

    code += "PIKG::S32 Retrieve#{$kernel_name}(const PIKG::S32 tag,\n"
    code += "                       const PIKG::S32 n_walk,\n"
    code += "                       const PIKG::S32 ni[],\n"
    code += "                       #{$force_name} *force[])\n"
    code += "{\n"
    code += "  int ni_tot = 0;\n"
    code += "  for(int k=0; k<n_walk; k++){\n"
    code += "    ni_tot += ni[k];\n"
    code += "  }\n"
    code += "  dev_force.d2h(ni_tot);\n"
    code += "\n"
    code += "  int n_cnt = 0;\n"
    code += "  for(int iw=0; iw<n_walk; iw++){\n"
    code += "    for(int i=0; i<ni[iw]; i++){\n"
    # accum force
    fvars.each{ |v|
      iotype = h[v][0]
      if iotype == "FORCE"
        type = h[v][1]
        modifier = h[v][3]
        if modifier == "local"
          fdpsname = v
        else
          fdpsname = h[v][2]
        end
        type_single = get_single_element_type(type)
        get_vector_elements(type).zip($accumhash[v]){ |dim,op|
          dest = "force[iw][i]." + fdpsname
          dest = Expression.new([:dot,dest,dim,type_single]) if type =~ /vec/
          src = "dev_force[n_cnt]." + fdpsname
          src = Expression.new([:dot,src,dim,type_single]) if type =~ /vec/
          code += "      " + Accumulate.new([dest,src,1,type_single,op]).convert_to_code("reference") + "\n"
        }
      end
    }

    code += "      n_cnt++;\n"
    code += "    }\n"
    code += "  }\n"
    code += "  return 0;\n"
    code += "}\n"

    code += "void clearTimeProfiler(){\n"
    code += "  prof.clearAll();\n"
    code += "}\n"


    File.open($output_file, mode = 'w'){ |f|
      f.write(code)
    }
  end
end

