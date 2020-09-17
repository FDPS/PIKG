class Load
  attr_accessor :dest, :src, :nelem, :type, :iotype
  def initialize(x)
    @dest,@src,@nelem,@type,@iotype = x
  end

  def convert_to_code(conversion_type)
    #abort "type of lval and rval is different" if @src.get_type != @dest.get_type
    ret = String.new

    tot, max_byte_size, is_uniform = count_class_member(@iotype)
    type_single = get_single_element_type(@type)
    size = get_single_data_size(type_single)
    lane_size = get_simd_width(conversion_type) / size
    nlane = lane_size / @nelem
    offset_gather = ((tot+max_byte_size-1)/max_byte_size * max_byte_size) / byte_count(type_single)
    case conversion_type
    when "reference"
      get_vector_elements(@type).each{ |dim|
        src = @src
        src = Expression.new([:dot,src,dim,get_single_element_type(@type)]) if dim != ""

        dest = @dest
        dest = Expression.new([:dot,dest,dim,get_single_element_type(@type)]) if dim != ""
        ret += Statement.new([dest,src,@type,nil]).convert_to_code(conversion_type)
      }

    when "AVX2"
      case nlane
      when lane_size
        get_vector_elements(@type).each{ |dim|
          suffix = get_type_suffix_avx2(type)
          src = @src
          src = Expression.new([:dot,src,dim,get_single_element_type(@type)]) if dim != ""
          dest = @dest
          dest = Expression.new([:dot,dest,dim,get_single_element_type(@type)]) if dim != ""

          ret += "#{@dest.convert_to_code(conversion_type)} = _mm256_set1_#{suffix}(#{src.convert_to_code("reference")});\n"
        }
      when 1
        ret += GatherLoad.new([dest,src,"0","#{offset_gather}",type]).convert_to_code(conversion_type)
      else
        abort "unsupported number of elemet (#{@nelem}) for Load" if nlane <= 0
        index = String.new
        if @iotype == "EPI" || @iotype == "FORCE"
          for j in 0...@nelem
            for i in 0...nlane
              index += "," if !(i == 0 && j == 0)
              index += "#{i*offset_gather}"
            end
          end
        elsif @iotype == "EPJ"
          for i in 0...nlane
            for j in 0...@nelem
              index += "," if !(i==0 && j==0)
              index += "#{i*offset_gather}"
            end
          end
        else
          abort
        end

        ret += GatherLoad.new([dest,src,nil,nil,type]).convert_to_code(conversion_type,index)
      end
    when "AVX-512"
      suffix = get_type_suffix_avx512(type)
      case nlane
      when lane_size
        get_vector_elements(@type).each{ |dim|
          src = @src
          src = Expression.new([:dot,src,dim,get_single_element_type(@type)]) if dim != ""
          dest = @dest
          dest = Expression.new([:dot,dest,dim,get_single_element_type(@type)]) if dim != ""

          ret += "#{@dest.convert_to_code(conversion_type)} = _mm512_set1_#{suffix}(#{src.convert_to_code("reference")});\n"
        }
      when 1
        ret += GatherLoad.new([dest,src,"0","#{offset_gather}",type]).convert_to_code(conversion_type)
      else
      end
    when "A64FX"
      suffix = get_type_suffix_a64fx(type)
      case nlane
      when lane_size
        get_vector_elements(@type).each{ |dim|
          src = @src
          src = Expression.new([:dot,src,dim,get_single_element_type(@type)]) if dim != ""
          dest = @dest
          dest = Expression.new([:dot,dest,dim,get_single_element_type(@type)]) if dim != ""

          ret += "#{@dest.convert_to_code(conversion_type)} = _svld1_#{suffix}(#{$current_predicate},#{src.convert_to_code("reference")});\n"
        }
      when 1
        ret += GatherLoad.new([dest,src,"0","#{offset_gather}",type]).convert_to_code(conversion_type)
      else
      end
    end
    ret
  end
end

class Accumulate
  attr_accessor :dest, :src, :nelem, :type, :op
  def initialize(x)
    @dest,@src,@nelem,@type,@op = x
  end
  def convert_to_code(conversion_type)
    ret = String.new

    tot, max_byte_size, is_uniform = count_class_member("FORCE")
    type_single = get_single_element_type(@type)
    size = get_single_data_size(type_single)
    lane_size = get_simd_width(conversion_type) / size
    nlane = lane_size / @nelem

    offset_scatter = ((tot+max_byte_size-1)/max_byte_size * max_byte_size) / byte_count(type_single)

    case conversion_type
    when "reference"
      get_vector_elements(@type).each{ |dim|
        dest = @dest
        dest = Expression.new([:dot,dest,dim,get_single_element_type(@type)]) if dim != ""
        src = @src
        src = Expression.new([:dot,src,dim,get_single_element_type(@type)]) if dim != ""
        if @op == "max" || @op == "min"
          ret += Statement.new([dest,FuncCall.new([@op,[src,dest],@type])]).covert_to_code(conversion_type)
        else
          ret += Statement.new([dest,src,@type,:plus]).convert_to_code(conversion_type) if @op == :plus || @op == :minus
          ret += Statement.new([dest,src,@type,:mult]).convert_to_code(conversion_type) if @op == :mult || @op == :div
        end
      }
    when "AVX2"
      case nlane
      when lane_size
        ret += "{\n"
        if @op == :plus || @op == :minus
          op = "add"
        elsif @op == :mult || @op == :div
          op = "mul"
        elsif @op == "max"
          op = "max"
        elsif @op == "min"
          op = "min"
        end
        src = @src.convert_to_code(conversion_type)
        suffix = get_type_suffix_avx2(@type)
        imm8 = "0b1111" if size == 64
        imm8 = "0xb1" if size == 32
        ret += "#{src} = _mm256_#{op}_#{suffix}(#{src},_mm256_shuffle_#{suffix}(#{src},#{src},#{imm8}));\n"
        ret += "#{src} = _mm256_#{op}_#{suffix}(#{src},_mm256_shuffle_#{suffix}(#{src},#{src},0xee));\n" if @type == "F32"
        ret += "#{src} = _mm256_#{op}_#{suffix}(#{src},_mm256_cast#{suffix}128_#{suffix}256(_mm256_extractf128_#{suffix}(#{src},1)));\n"
        dest_conv = "#{@dest.convert_to_code(conversion_type)}[0]"
        src_conv = "#{src}[0]"
        if op == "max" || op == "min"
          ret += "#{dest_conv} = #{op}(#{dest_conv},#{src_conv});"
        else
          case op
          when "add"
            ret += NonSimdState.new([dest_conv,NonSimdExp.new([:plus,dest_conv,src_conv,@type]),@type,nil]).convert_to_code(conversion_type) + "\n"
          when "mul"
            ret += NonSimdState.new([dest_conv,NonSimdExp.new([:mult,dest_conv,src_conv,@type]),@type,nil]).convert_to_code(conversion_type) + "\n"
          end
        end
        ret += "}\n"
      when 1
        tmp = "__fkg_tmp_accum"
        ret += "{\n"
        ret += Declaration.new([type,tmp]).convert_to_code(conversion_type)
        ret += GatherLoad.new([tmp,dest,"0","#{offset_scatter}",type]).convert_to_code(conversion_type) + "\n"
        if @op == "max" || @op == "min"
          rexp = FuncCall.new([@op,[tmp,@src],type])
        elsif @op == :plus || @op == :minus
          rexp = Expression.new([:plus,tmp,@src,type])
        elsif @op == :mult || @op == :div
          rexp = Expression.new([:mult,tmp,@src,type])
        end
        ret += Statement.new([tmp,rexp]).convert_to_code(conversion_type) + "\n"
        ret += ScatterStore.new([dest,tmp,"0","#{offset_scatter}",type]).convert_to_code(conversion_type)
        ret += "}\n"
      else
        abort "Accumulate of nlane #{nlane} for AVX2 is not supported"
      end
    when "AVX-512"
      if @op == :plus || @op == :minus
        op = "add"
      elsif @op == :mult || @op == :div
        op = "mul"
      elsif @op == "max"
        op = "max"
      elsif @op == "min"
        op = "min"
      end
      suffix = get_type_suffix_avx512(@type)
      src = @src.convert_to_code(conversion_type)
      case nlane
      when lane_size # 16 for FP32, 8 for FP64
        ret += "#{@dest.convert_to_code(conversion_type)}[0] += _mm512_reduce_#{op}_#{suffix}(#{src});\n"
      when 1
        tmp = "__fkg_tmp_accum"
        ret += "{\n"
        ret += Declaration.new([type,tmp]).convert_to_code(conversion_type)
        ret += GatherLoad.new([tmp,dest,"0","#{offset_scatter}",type]).convert_to_code(conversion_type) + "\n"
        if @op == "max" || @op == "min"
          rexp = FuncCall.new([@op,[tmp,@src],type])
        elsif @op == :plus || @op == :minus
          rexp = Expression.new([:plus,tmp,@src,type])
        elsif @op == :mult || @op == :div
          rexp = Expression.new([:mult,tmp,@src,type])
        end
        ret += Statement.new([tmp,rexp]).convert_to_code(conversion_type) + "\n"
        ret += ScatterStore.new([dest,tmp,"0","#{offset_scatter}",type]).convert_to_code(conversion_type)
        ret += "}\n"
      when "A64FX"
        if @op == :plus || @op == :minus
          op = "add"
        elsif @op == :mult || @op == :div
          op = "mul"
        elsif @op == "max"
          op = "max"
        elsif @op == "min"
          op = "min"
        end
        suffix = get_type_suffix_avx512(@type)
        src = @src.convert_to_code(conversion_type)
        case nlane
        when lane_size # 16 for FP32, 8 for FP64
          ret += "#{@dest.convert_to_code(conversion_type)}[0] += _mm512_reduce_#{op}_#{suffix}(#{src});\n"
        when 1
          tmp = "__fkg_tmp_accum"
          ret += "{\n"
          ret += Declaration.new([type,tmp]).convert_to_code(conversion_type)
          ret += GatherLoad.new([tmp,dest,"0","#{offset_scatter}",type]).convert_to_code(conversion_type) + "\n"
          if @op == "max" || @op == "min"
            rexp = FuncCall.new([@op,[tmp,@src],type])
          elsif @op == :plus || @op == :minus
            rexp = Expression.new([:plus,tmp,@src,type])
          elsif @op == :mult || @op == :div
            rexp = Expression.new([:mult,tmp,@src,type])
          end
          ret += Statement.new([tmp,rexp]).convert_to_code(conversion_type) + "\n"
          ret += ScatterStore.new([dest,tmp,"0","#{offset_scatter}",type]).convert_to_code(conversion_type)
          ret += "}\n"
        else
          abort "unsupported number of lane per element for A64FX"
        end
      else
        abort "unnsupported conversion_type #{coversion_type} for Accumulate"
      end
    end
    ret
  end
end

  class TailJLoop
    attr_accessor :loop,:ninj,:fvars
    def initialize(x)
      @loop, @ninj, @fvars = x
    end

    def convert_to_code(conversion_type)
      ret = String.new
      ret += "if(j<nj){ // tail j loop\n"
      tmp_hash = Hash.new
      fvars.each{ |name,type|
        size_single = get_single_data_size(type)
        n_split =  size_single / $min_element_size
        for i in 0...n_split
          tmp = add_new_tmpvar(type)
          ret += Declaration.new([type,tmp]).convert_to_code(conversion_type) + "\n"
          type_single = get_single_element_type(type)
          get_vector_elements(type).each{ |dim|
            src = name
            src += "_#{i}" if n_split > 1
            tmp_hash[src] = tmp if tmp_hash[src] == nil
            src = Expression.new([:dot,src,dim,type_single]) if dim != ""
            dest = tmp
            dest = Expression.new([:dot,dest,dim,type_single]) if dim != ""
            ret += Statement.new([dest,src,type_single,nil]).convert_to_code(conversion_type) + "\n"
          }
        end
      }
      ret += @loop.convert_to_code(conversion_type)
      case conversion_type
      when "refference"
        abort "reference mode does not reach TailJLoop"
      when "AVX2"
        fvars.each{ |name,type|
          size_single = get_single_data_size(type)
          n_split =  size_single / $min_element_size
          type_single = get_single_element_type(type)
          case @ninj[1]
          when 1
            abort "Not need tail loop of j"
          when 2
            imm8 = "0b00000011" if type_single =~ /64/
            imm8 = "0b00001111" if type_single =~ /32/
          when 4
            imm8 = "0b00000001" if type_single =~ /64/
            imm8 = "0b00000011" if type_single =~ /32/
          when 8
            abort "j_parallel of F64 must be <= 4" if type_single == /64/
            imm8 = "0b00000001" if type_single == "F32"
          end
          for i in 0...n_split
            get_vector_elements(type).each{ |dim|
              dest = name
              dest += "_#{i}" if n_split > 1
              src = tmp_hash[dest]
              dest = Expression.new([:dot,dest,dim,type_single]) if dim != ""
              src = Expression.new([:dot,src,dim,type_single]) if dim != ""
              suffix = get_type_suffix_avx2(type_single)
              ret += "#{dest.convert_to_code(conversion_type)} = _mm256_blend_#{suffix}(#{src.convert_to_code(conversion_type)},#{dest.convert_to_code(conversion_type)},#{imm8});\n"
            }
          end
        }
      when "AVX-512"
        fvars.each{ |name,type|
          size_single = get_single_data_size(type)
          n_split =  size_single / $min_element_size
          type_single = get_single_element_type(type)
          case @ninj[1]
          when 1
            abort "Not need tail loop of j"
          when 2
            imm8 = "0b00001111" if type_single =~ /64/
            imm8 = "0b11111111" if type_single =~ /32/
          when 4
            imm8 = "0b00000011" if type_single =~ /64/
            imm8 = "0b00001111" if type_single =~ /32/
          when 8
            imm8 = "0b00000001" if type_single =~ /F64/
            imm8 = "0b00000011" if type_single =~ /F32/
          when 16
            abort "j_parallel of F64 must be <= 8" if type_Single == /64/
            imm8 = "0b00000001" if type_single =~ "F32"
          end
          for i in 0...n_split
            get_vector_elements(type).each{ |dim|
              dest = name
              dest += "_#{i}" if n_split > 1
              src = tmp_hash[dest]
              dest = Expression.new([:dot,dest,dim,type_single]) if dim != ""
              src = Expression.new([:dot,src,dim,type_single]) if dim != ""
              suffix = get_type_suffix_avx2(type_single)
              ret += "#{dest.convert_to_code(conversion_type)} = _mm512_mask_blend_#{suffix}(_cvtu32_mask#{get_simd_width(conversion_type)/size_single}(#{imm8}),#{src.convert_to_code(conversion_type)},#{dest.convert_to_code(conversion_type)});\n"
            }
          end
        }
      when "A64FX"
        abort "conversion_type A64FX does not need tail loop"
      end
      ret += "} // if of j tail loop\n"
      ret
    end
  end

  class SplitPredicate
    attr_accessor :nsplit
    def initialize(x)
      @nsplit = x
    end
    def convert_to_code(conversion_type)
      ret = String.new
      case conversion_type
      when "reference"
      # do nothing
      when "A64FX"
        orig = $current_predicate
        for i in 0...@nsplit
          name = orig + "_#{i}"
          ret += "svbool_t #{name};\n"
        end
      when "AVX2"
        orig = $current_predicate
        for i in 0...@nsplit
          name = orig + "_#{i}"
          ret += "__m256 #{name};\n"
        end
        case @nsplit
        when 2
          ret += "#{orig}_0 = _mm256_castpd_ps(_mm256_cvtps_pd(_mm256_castps256_ps128(#{orig})));\n"
          ret += "#{orig}_1 = _mm256_castpd_ps(_mm256_cvtps_pd(_mm256_extractf128_ps(#{orig},1)));\n"
        when 4
          abort "unsupported n_split for SplitPredicate of AVX2"
        end
      when "AVX-512"
        orig = $current_predicate
        for i in 0...@nsplit
          name = orig + "_#{i}"
          ret += "__mmask#{get_simd_width(conversion_type)/$min_element_size/@nsplit} #{name};\n"
        end
        case @nsplit
        when 2
          ret += "#{orig}_0 = _cvtu32_mask8(_cvtmask16_u32(#{orig}));\n"
          ret += "#{orig}_1 = _cvtu32_mask8(_cvtmask16_u32(#{orig})>>8);\n"
        when 4
          abort "unsupported n_split for SplitPredicate of AVX-512"
        end
      end
      ret
    end
  end

  class Kernelprogram
    def generate_optimized_code_multi_prec(conversion_type,output = $output_file)
      code = "#include<pikg_vector.hpp>\n"
      code +="#include<cmath>\n"
      code +="#include<limits>\n"
      code +="#include<chrono>\n"
      code += $additional_text if $additional_text != nil
      code += "\n"
      case conversion_type
      when /A64FX/
        code += "#include <arm_sve.h>\n"
      when /AVX2/
        code += "#include <pikg_avx2.hpp>\n"
      when /AVX-512/
        code += "#include <pikg_avx512.hpp>\n"
      end

      if $c_interface_impl
        struct_list = ["EPI"]
        struct_list.push("EPJ") if $epi_name != $epj_name
        struct_list.push("FORCE") if $epi_name != $force_name && $epj_name != $force_name
        struct_list.zip([$epi_name,$epj_name,$force_name]){ |c,n|
          next if n.index("PS::")
          code +="struct #{n}{\n"
          $varhash.each{|v|
            iotype = get_iotype_from_hash(v);
            modifier = get_modifier_from_hash(v)
            if iotype == c && modifier == nil
              type = get_type_from_hash(v)
              fdpsname = get_fdpsname_from_hash(v)
              code += "PIKG::#{type} #{fdpsname};\n"
            end
          }
          code += "};\n"
        }
      end


      fvars = generate_force_related_map(@statements)
      # calc lane_size
      $min_element_size = 64
      $max_element_size = 16
      fvars.each{ |v|
        element_size = get_single_data_size($varhash[v][1])
        $min_element_size = element_size if element_size < $min_element_size
        $max_element_size = element_size if element_size > $max_element_size
      }

      ninj_set = []
      ni = get_simd_width(conversion_type) / $max_element_size  #lane_size
      ni = 1 if ni <= 0
      nj = 1
      ninj_set.push([ni,nj])
      ninj_set.push([nj,ni]) if ni > 1
      #while ni > 0
      #  ninj_set += [[ni,nj]]
      #  ni /= 2
      #  nj *= 2
      #  #abort "ni(#{ni}) * nj(#{nj}) != lane_size(#{lane_size})" if ni*nj != lane_size
      #end
      #p ninj_set

      code += kernel_class_def(conversion_type)
      code += "if(kernel_select>=0) kernel_id = kernel_select;\n"
      code += "if(kernel_id == 0){\n"
      code += "#{$force_name}* force_tmp = new #{$force_name}[ni];\n"
      code += "std::chrono::system_clock::time_point  start, end;\n"
      code += "double min_time = std::numeric_limits<double>::max();\n"
      ninj_set.each_with_index{ |ninj,index|
        code += "{ // test Kernel_I#{ninj[0]}_J#{ninj[1]}\n"
        code += "for(int i=0;i<ni;i++) force_tmp[i] = force[i];\n"
        code += "start = std::chrono::system_clock::now();\n"
        code += "Kernel_I#{ninj[0]}_J#{ninj[1]}(epi,ni,epj,nj,force_tmp);\n"
        code += "end = std::chrono::system_clock::now();\n"
        code += "double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();\n"
        #code += "std::cerr << \"kerel #{index+1}: \" << elapsed << \" ns\" << std::endl;\n"
        code += "if(min_time > elapsed){\n"
        code += "min_time = elapsed;\n"
        code += "kernel_id = #{index+1};\n"
        code += "}\n"
        code += "}\n"
      }
      code += "delete[] force_tmp;\n"
      #code += "std::cerr << \"kernel \" << kernel_id << \" is selected\" << std::endl;\n"
      code += "} // if(kernel_id == 0)\n"

      ninj_set.each_with_index{ |ninj,index|
        code += "if(kernel_id == #{index+1}) Kernel_I#{ninj[0]}_J#{ninj[1]}(epi,ni,epj,nj,force);\n"
      }
      code += "} // operator() definition \n"

      ninj_set.each{|ninj|
        $varhash.each{|v|
          iotype   = v[1][0]
          v[1][0] = nil if iotype == "declared"
        }
        code += "void Kernel_I#{ninj[0]}_J#{ninj[1]}(const #{$epi_name}* __restrict__ epi,const PIKG::S32 ni,const #{$epj_name}* __restrict__ epj,const PIKG::S32 nj,#{$force_name}* __restrict__ force){\n"
        $pg_count = 0
        code += kernel_body_multi_prec(ninj,conversion_type)
        code += "} // Kernel_I#{ninj[0]}_J#{ninj[1]} definition \n"
      }

      code += reserved_func_def(conversion_type)
      code += "};// kernel functor definition \n"

      if $c_interface_impl
        code += "#{$kernel_name} __pikg_#{$interface_name};\n"

        code += "extern \"C\"{\n"
        code += "void #{$initializer_name}("
        count = 0
        $varhash.each{|v|
          iotype = v[1][0]
          if iotype == "MEMBER"
            code += "," if count > 0
            name = v[0]
            type = v[1][1]
            code += "PIKG::" + type + " " + name +"_"
            count = count + 1
          end
        }
        code += "){\n"
        code += "__pikg_#{$interface_name}.initialize("
        count = 0
        $varhash.each{|v|
          iotype = v[1][0]
          if iotype == "MEMBER"
            code += "," if count > 0
            name = v[0]
            type = v[1][1]
            code += name +"_"
            count = count + 1
          end
        }
        code += ");\n"
        code += "} // intialize_#{$interface_name}\n"

        code += "void #{$interface_name}(const #{$epi_name}* epi, const PIKG::S32 ni, const #{$epj_name}* epj,const PIKG::S32 nj, #{$force_name}* force){\n"
        code += "__pikg_#{$interface_name}(epi,ni,epj,nj,force);\n"
        code += "}\n"
        code += "} // extern \"C\"\n"
      end
      
      File.open(output, mode = 'w'){ |f|
        f.write(code)
      }
    end

    def generate_loop_begin_multi_prec(conversion_type,lane_size,ij,start = nil)
      n = "n"+ij
      case conversion_type
      when "reference"
        ret = Loop.new([ij,start,n,1,[]])
      when  /A64FX/
        ret = Loop.new([ij,start,"((#{n}+#{lane_size-1})/#{lane_size})*#{lane_size}","#{lane_size}",["#{$current_predicate}=svwhilelt_b32_s32(#{ij},#{n});"]])
      when  /AVX/
        ret = Loop.new([ij,start,"(#{n}/#{lane_size})*#{lane_size}","#{lane_size}",[]])
      end
      ret
    end

    class Fusion
      attr_accessor :ops,:from, :to
      def initialize(x)
        @ops, @from, @to = x
        abort "type of Fusion must be specified" if @from == nil || @to == nil
      end
      def convert_to_code(conversion_type)
        ret = String.new
        n = @ops.length
        case conversion_type
        when "reference"
        # do nothing
        when "AVX2"
          case n
          when 1
            abort "Fusion does not happen for single operand"
          when 2
            if @from == "F64" && @to == "F32"
              ret = "_mm256_permute2f128_ps("
              for i in 0...n
                ret += "," if i>0
                ret += "_mm256_castps128_ps256(_mm256_cvtpd_ps(#{ops[i]}))"
              end
              ret += ",0x20)"
            else
              abort "Fusion for AVX2 under construction"
            end
          else
            abort "Fusion of #{n} does not happen for AVX2"
          end

        when "AVX-512"
          case n
          when 1
            abort "Fusion does not happen for single operand"
          when 2
            if @from == "F64" && @to == "F32"
              ret = "_mm512_insertf32x8(_mm512_castps256_ps512(_mm512_cvtpd_ps(#{ops[0]})),_mm512_cvtpd_ps(#{ops[1]}),1)"
            end
          else
            abort "Fusion #{n} for AVX2 under construction"
          end
        else
          abort "unsupported conversion_type #{coversion_type} for Fusion"
        end
        ret
      end
    end

    class Fission
      attr_accessor :op,:from, :to, :index
      def initialize(x)
        @op, @from, @to, @index = x
        abort "type of Fision must be specified" if @from == nil || @to == nil
        abort "index must be specified" if @index == nil
      end
      def convert_to_code(conversion_type)
        ret = String.new
        case conversion_type
        when "reference"
        # do nothing
        when "AVX2"
          if @from == "F32" && @to == "F64"
            ret = @op.convert_to_code(conversion_type)
            ret = "_mm256_castps256_ps128(#{ret})"  if @index == 0
            ret = "_mm256_extractf128_ps(#{ret},1)" if @index == 1
            ret = "_mm256_cvtps_pd(#{ret})"
          else
            abort "under costruction for AVX2 : Fission"
          end
        when "AVX-512"
          if @from == "F32" && @to == "F64"
            ret = @op.convert_to_code(conversion_type)
            ret = "_mm512_castps512_ps256(#{ret})"  if @index == 0
            ret = "_mm512_extractf32x8_ps(#{ret},1)" if @index == 1
            ret = "_mm512_cvtps_pd(#{ret})"
          else
            abort "under costruction for AVX2 : Fission"
          end
        else
          abort "unsupported conversion_type #{conversion_type} for Fission"
        end
        ret
      end
    end

    def split_downcast(ss)
      reet = Array.new

      downcasting = false
      to = nil
      from = nil
      nline = nil
      ss.reverse_each{ |s|
        if s.expression.class == FuncCall && s.expression.name =~ /to\_(f|s|u)(64|32|16)/
          to = s.name.get_type
          from = s.expression.ops[0].get_type
          downcasting = true
          #ret += ["//downcasting from #{from} to #{to}"]

          nline = get_single_data_size(from) / get_single_data_size(to)
          p "nline:",nline
          ops = Array.new
          tmps = Array.new
          for i in 0...nline
            orig = get_name(s)
            replace = "#{orig}_#{i}"
            tmp = Statement.new([s.name,s.expression.ops[0],s.type,s.op])
            tmp.replace_name(orig,replace)

            tmps += [tmp]
            ops += [replace]
          end
          ret += [Statement.new([get_name(s),Fusion.new([ops,from,to]),to,:nil])]
          ret += [s.declare_temporal_var]
          ret += tmps.reverse
        elsif downcasting
          for i in 0...nline
            ret += [s]
          end
        elsif s.class == ConditioalBranch
          tmp_cbb = ConditionalBranch.new([[],[]])
          tmp_cbb.conditions = s.conditions
          s.bodies.each{ |bss|
            tmp_cbb.bodies.push(split_downcast(bss))
          }
          ret += [tmp_cbb]
        else
          ss += [s]
        end
      }
      ret = ss.reverse
      ret
    end

    def load_local_var(name,type,offset = 0)
      get_vector_elements(type).each{|dim|
        dst = Expression.new([:dot,name,dim,tmp_type])
        src = PointerOf.new([tmp_type,Expression.new([:array,"#{name}_tmp_"+dim,"i+#{offset}",])])
        ret += [LoadState.new([dst,src,tmp_type])]
      }
      ret
    end

    def load_vars(io,fvars,nelem,conversion_type,h=$varhash)
      ret = Array.new
      fvars.each{ |v|
        iotype = h[v][0]
        if iotype == io
          type = h[v][1]
          fdpsname = h[v][2]
          modifier = h[v][3]
          n = get_single_data_size(type) / $min_element_size
          n = 1 if conversion_type == "reference"
          for i in 0...n
            name = v
            name += "_#{i}" if n > 1
            offset = i*nelem #get_num_elem(type,conversion_type)
            ret += [Declaration.new([type,name])]
            if modifier == "local"
              ret += load_local_var(name,type,offset)
            else
              tot, max_byte_size, is_uniform = count_class_member(iotype)
              type_single = get_single_element_type(type)
              offset_gather = ((tot+max_byte_size-1)/max_byte_size * max_byte_size) / byte_count(type_single)
              get_vector_elements(type).each{ |dim|
                index = "i+#{offset}" if io == "EPI" || io == "FORCE"
                index = "j+#{offset}" if io == "EPJ"
                src = Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i+#{offset}"]),fdpsname,type_single])
                src = Expression.new([:dot,src,dim,type_single]) if type =~ /vec/
                src = PointerOf.new([type,src])
                dest = name
                dest = Expression.new([:dot,name,dim,type_single]) if type =~ /vec/
                ret += [Load.new([dest,src,nelem,type_single,iotype])]
              }
            end
            h[name] = [iotype,type,fdpsname,"alias"] if h[name] == nil
          end
        end
      }
      ret
    end

    def load_jvars(fvars,nelem,conversion_type,h=$varhash)
      ret = Array.new
      fvars.each{ |v|
        iotype = h[v][0]
        if iotype == "EPJ"
          type = h[v][1]
          fdpsname = h[v][2]
          modifier = h[v][3]

          name = v
          ret += [Declaration.new([type,name])]
          if modifier == "local"
            ret += load_local_var(name,type,offset)
          else
            tot, max_byte_size, is_uniform = count_class_member(iotype)
            type_single = get_single_element_type(type)
            offset_gather = ((tot+max_byte_size-1)/max_byte_size * max_byte_size) / byte_count(type_single)
            get_vector_elements(type).each{ |dim|
              index = "j"
              src = Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),index]),fdpsname,type_single])
              src = Expression.new([:dot,src,dim,type_single]) if type =~ /vec/
              src = PointerOf.new([type,src])
              dest = name
              dest = Expression.new([:dot,name,dim,type_single]) if type =~ /vec/
              ret += [Load.new([dest,src,nelem,type_single,iotype])]
            }
          end
        end
      }
      ret
    end

    def init_force(fvars,accum_hash,conversion_type,h=$varhash)
      ret = Array.new
      fvars.each{ |v|
        iotype = h[v][0]
        if iotype == "FORCE"
          type = h[v][1]
          fdpsname = h[v][2]
          type_single = get_single_element_type(type)
          n = get_single_data_size(type) / $min_element_size
          n = 1 if conversion_type == "reference"
          for i in 0...n
            name = v
            name += "_#{i}" if n > 1
            ret += [Declaration.new([type,name])]
            get_vector_elements(type).each_with_index{ |dim,j|
              dest = name
              dest = Expression.new([:dot,dest,dim,type_single]) if dim != ""
              op = accum_hash[v][j]
              abort "accum_hash == nil" if accum_hash[v][j] == nil
              ret += [Duplicate.new([dest,get_initial_value(op,type_single),type_single])]
            }
            h[name] = [iotype,type,fdpsname,"alias"] if h[name] == nil
          end
        end
      }
      ret
    end
    
    def store_vars(accum_hash,fvars,nelem,conversion_type,h=$varhash)
      ret = Array.new
      fvars.each{ |v|
        iotype = h[v][0]
        if iotype == "FORCE"
          type = h[v][1]
          fdpsname = h[v][2]
          modifier = h[v][3]
          n = get_single_data_size(type) / $min_element_size
          n = 1 if conversion_type == "reference"
          for i in 0...n
            name = v
            name += "_#{i}" if n > 1
            offset = i*nelem #get_num_elem(type,conversion_type)

            type_single = get_single_element_type(type)
            get_vector_elements(type).zip(accum_hash[v]){ |dim,op|
              dest = Expression.new([:dot,Expression.new([:array,get_iotype_array(iotype),"i+#{offset}"]),fdpsname,type_single])
              dest = Expression.new([:dot,dest,dim,type_single]) if type =~ /vec/
              dest = PointerOf.new([type,dest])
              src = name
              src = Expression.new([:dot,name,dim,type_single]) if type =~ /vec/
              ret += [Accumulate.new([dest,src,nelem,type_single,op])]
            }
          end
        end
      }
      ret
    end

    def generate_jloop_body(ss,fvars,split_vars,conversion_type,h=$varhash)
      ret = Array.new
      ss.each{ |s|
        if conversion_type == "reference"
          if isStatement(s)
            ret += s.declare_temporal_var 
            ret.push(s)
          elsif s.class == ConditionalBranch
            tmp_cbb = ConditionalBranch.new([[],[]])
            tmp_cbb.conditions = s.conditions
            s.bodies.each{ |bss|
              tmp_cbb.bodies.push(generate_jloop_body(bss,fvars,split_vars,conversion_type,h))
            }
            ret += [tmp_cbb]
          else
            ret += s.declare_temporal_var
            ret += [s]
          end
        else
          if isStatement(s)
            lexp = get_name(s)
            list = [lexp] + s.expression.get_related_variable
            doSplit = false
            list.each{ |v| doSplit = true if split_vars.index(v)}
            if s.expression.class == FuncCall && s.expression.name =~ /to_(f|s|u)(64|32|16)/
              #p split_vars
              type_to = s.get_type
              type_from = s.expression.ops[0].get_type
              size_to = get_single_data_size(type_to)
              size_from= get_single_data_size(type_from)
              n_to   = size_to   / $min_element_size
              n_from = size_from / $min_element_size
              if size_to == size_from
                ret += [Statement.new([s.name,s.expression.ops[0],s.type,s.op])]
              elsif size_to < size_from
                #p "# of fusion var:",n_to,n_from
                ops = Array.new
                for i in 0...n_from
                  name = lexp + "_#{i}"
                  tmp = Statement.new([s.name.dup,s.expression.ops[0].dup,type_from,s.op])
                  tmp.replace_name(lexp,name)
                  #p s
                  split_vars.each{ |orig|
                    #p "split:",orig
                    #p tmp.expression
                    tmp.expression = tmp.expression.replace_recursive(orig,orig + "_#{i}")
                  }
                  ret += [Declaration.new([type_from,name])]
                  ret += [tmp]
                  ops += [name]
                end
                ret += s.declare_temporal_var
                ret += [Statement.new([lexp,Fusion.new([ops,type_from,type_to]),type_to])]
              elsif size_to > size_from
                # fission
                rexp = s.expression.ops[0]
                for i in 0...n_to
                  op = rexp
                  op += "_#{i/n_from}" if n_from>1
                  name = lexp + "_#{i}"
                  tmp = Statement.new([name,Fission.new([op,type_from,type_to,i%(n_to/n_from)]),type_to])

                  $varhash[name] = [nil,type_to,nil] if $varhash[name] == nil
                  ret += tmp.declare_temporal_var
                  ret += [tmp]
                end
                split_vars += [lexp]
              else
                ret += [s]
              end
            elsif doSplit
              type = s.type
              n = get_single_data_size(type) / $min_element_size
              n = 1 if conversion_type == "reference"
              if n > 1
                for i in 0...n
                  name = lexp + "_#{i}"
                  tmp = Statement.new([s.name.dup,s.expression.dup,s.type,s.op])
                  tmp.replace_name(lexp,name)
                  split_vars.each{ |orig|
                    tmp.expression = tmp.expression.replace_recursive(orig,orig + "_#{i}")
                  }
                  tmp.expression.split_index = i if tmp.expression.class == Merge

                  $varhash[name] = [nil,@type,nil,nil] if $varhash[name] == nil
                  ret += tmp.declare_temporal_var($varhash)
                  ret += [tmp]

                end
              end

              split_vars += [lexp]
            else
              ret += s.declare_temporal_var
              ret += [s]
            end
          elsif s.class == ConditionalBranch
            tmp_cbb = ConditionalBranch.new([[],[]])
            tmp_cbb.conditions = s.conditions
            tmp_predicate = Array.new
            s.bodies.zip(s.conditions){ |bss,cond|
              tmp = Array.new
              tmp_nsplit = Array.new
              bss.each{ |bs|
                if fvars.index(get_name(bs)) != nil
                  n = get_single_data_size(bs.type) / $min_element_size
                  if n > 1
                    tmp_nsplit.push(n)
                  end
                end
              }
              tmp_nsplit.sort.uniq.each{ |n|
                tmp += [SplitPredicate.new(n)]
              }
              tmp_predicate.push(tmp)
            }
            
            s.bodies.zip(tmp_predicate){ |bss,pss|
              tmp_cbb.bodies.push(pss + generate_jloop_body(bss,fvars,split_vars,conversion_type,h))
            }
            ret += [tmp_cbb]
          else
            ret += s.declare_temporal_var
            ret += [s]
          end
        end
      }
      ret
    end


    def kernel_body_multi_prec(ninj,conversion_type,istart=0,h=$varhash)
      #accum_hash = generate_accum_hash(@statements,h)

      code = String.new

      # declare "local" variables
      h.each{ |v|
        modifier = v[1][3]
        if modifier == "local"
          name = v[0]
          iotype = v[1][0]
          type   = v[1][1]
          array_size = ["ni","nj"][["EPI","EPJ"].index(iotype)]
          if type =~ /vec/
            type = type.delete("vec")
            ret += [Declaration.new([type,Expression.new([:array,"#{name}_tmp_x",array_size,type])])]
            ret += [Declaration.new([type,Expression.new([:array,"#{name}_tmp_y",array_size,type])])]
            ret += [Declaration.new([type,Expression.new([:array,"#{name}_tmp_z",array_size,type])])]
          else
            ret += [Declaration.new([type,Expression.new([:array,"#{name}_tmp",array_size,type])])]
          end
        end
      }

      # calc or copy local variables from EPI, EPJ, or FORCE
      @statements.each{ |s|
        name = get_name(s) if s.class == Statement
        if name != nil && h[name][3] == "local"
          tail = get_tail(s.name)
          iotype = h[name][0]
          type   = h[name][1]
          exp = s.expression
          index = get_io_index(iotype)
          loop_tmp = Loop.new([index,"0","n#{index}",1,[]])
          if tail != nil
            new_name = Expression.new([:array,"#{name}_tmp_#{tail}",index,type])
          else
            new_name = Expression.new([:array,"#{name}_tmp",index,type]) #"#{name}_tmp[i]"
          end
          new_exp = exp.replace_fdpsname_recursive(h)
          loop_tmp.statements += [Statement.new([new_name,new_exp])]
          ret += [loop_tmp]
        end
      }

      # declare and load TABLE variable
      tmp = []
      @statements.each{ |s|
        if s.class == TableDecl
          ret.push(s)
          tmp.push(s)
        end
      }
      tmp.each{ |s|
        @statements.delete(s)
      }

      fvars = generate_force_related_map(@statements)
      lane_size = get_simd_width(conversion_type) / $min_element_size
      lane_size = 1 if lane_size == 0

      split_vars = Array.new
      ["EPI","FORCE"].each{ |io|
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
      #p "split_vars:"
      #p split_vars
      iloop = generate_loop_begin_multi_prec(conversion_type,ninj[0]*$max_element_size/$min_element_size,"i",istart)
      #iloop.statements += load_ivars(conversion_type,lane_size,lane_size/ninj[0],fvars)

      # load EPI and FORCE variable
      iloop.statements += load_vars("EPI",fvars,ninj[0],conversion_type)
      #iloop.statements += load_vars("FORCE",fvars,ninj[0],conversion_type)
      iloop.statements += init_force(fvars,$accumhash,conversion_type)

      jloop = generate_loop_begin_multi_prec(conversion_type,ninj[1],"j",0)
      jloop.statements += load_jvars(fvars,ninj[1],conversion_type)
      jloop.statements += generate_jloop_body(@statements,fvars,split_vars,conversion_type)
      iloop.statements += [jloop]
      if ninj[1] > 1
        h.each{ |v|
          iotype = get_iotype_from_hash(v)
          v[1][0] = nil if iotype == "declared"
        }
        tmpvar = Array.new
        fvars.each{ |v|
          iotype = h[v][0]
          if iotype == "FORCE"
            type = h[v][1]
            tmpvar.push([v,type])
          end
        }
        jloop_tail = generate_loop_begin_multi_prec(conversion_type,1,"j",nil)
        jloop_tail.statements += load_jvars(fvars,1,conversion_type)
        jloop_tail.statements += generate_jloop_body(@statements,fvars,split_vars,conversion_type,nil)
        jloop_tail = TailJLoop.new([jloop_tail,ninj,tmpvar])
        iloop.statements += [jloop_tail]
      end
      iloop.statements += store_vars($accumhash,fvars,ninj[0],conversion_type)

      ss = Array.new
      ss += [iloop]

      code += NonSimdDecl.new(["S32","i"]).convert_to_code(conversion_type) if istart == 0
      code += NonSimdDecl.new(["S32","j"]).convert_to_code(conversion_type)
      ss.each{ |s|
        code += s.convert_to_code(conversion_type)
      }
      if conversion_type == "AVX2" || conversion_type == "AVX-512"
        h.each{|v|
          iotype   = v[1][0]
          v[1][0] = nil if iotype == "declared"
        }
        code += "{ // tail loop of reference \n"
        code += kernel_body("reference",nil,h)
        code += "} // end loop of reference \n"
      end
      #abort
      code
    end
  end
