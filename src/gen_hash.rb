$type = /(const)?\ +((PS::|PIKG::)?((F|S)(64|32|16)(vec(2|3|4)?)?)|(double|float((64|32|16)_t)?|(unsigned\ +)?(long\ +)+?int((64|32|16)_t)?))/
$ident=/[a-zA-Z_][a-zA-Z_0-9]*/

class Kernelprogram
  def process_iodecl(h = $varhash)
    @iodeclarations.each{|x|
      h[x.name] = [x.iotype, x.type, x.fdpsname, x.modifier]
    }
  end

  def process_funcdecl(h = $funchash)
    @functions.each{|x|
      decl = x.decl
      stmt = x.statements
      ret  = x.retval
      h[decl.name] = x
    }
  end

  def fusion_iotag
    ["EPI","EPJ","FORCE"].each{ |iotag|
      @statements.each{|s|
        s.fusion_iotag(iotag)
      }
    }
  end

  def generate_hash_from_cpp(filename,iotype,class_name,h = $varhash)
    code = String.new
    File.open(filename){ |f|
      f.each_line{ |line|
        next if line[0] == '#'
        next if line =~ /^\/\//
        code += line.chomp
      }
    }
    code = code.gsub(";","\n").gsub("{","{\n").gsub("}","}\n").gsub("public:",'').gsub("private:",'')
    #warn code

    #warn "test result:\n"
    nest_level = 0
    ignore = 0
    base_level = 1
    read = false
    code.each_line{ |line|
      next if line =~ /using/
      base_level += 1 if line =~ /namespace/
      if line =~ /(class|struct)\ +#{class_name}/
        read = true
      end
      if read
        nest_level += 1 if line =~ /\{/
        if nest_level == base_level
          if line =~ /^\s+#{$type}\s+#{$ident}(\s*,\s*#{$ident})*\s*\n/
            tmp = line.split(/(\s|,)/).select{ |s| s=~/(#{$type}|#{$ident})/}
            type = String.new
            type = tmp.shift
            if type == "const"
              next if iotype == "FORCE"
              type = tmp.shift
            end
            type = type.gsub("PS::","").gsub("PIKG::","")
            vars = tmp
            vars.each{ |v|
              h[iotype+"."+v] = [iotype,type,v,nil]
            }
          elsif line =~ /#{$ident}\ +#{$ident}(\ *,\ *#{$ident})*\ *\n/
            warn "undefined type: #{line}"
          end
        end
        nest_level -= 1 if line =~ /\}/
        read = false if nest_level == 0
      end
    }
  end

  def append_iotag(statement,iotag,h=$varhash)
    rexps = statement.expression.get_related_variable
    rexps.each{ |rexp|
      type = nil
      h.each{ |v|
        iotype = get_iotype_from_hash(v)
        fdpsname = get_fdpsname_from_hash(v)
        if iotype == iotag && fdpsname == rexp
          type = get_type_from_hash(v)
        end
      }
      statement.expression.replace_recursive(rexp,iotag+"_"+rexp) if type != nil
    }
  end
  
  def generate_alias(h=$varhash)
    fusion_iotag
    process_iodecl_alias(@iodeclarations,h)
    new_h = Hash.new
    h.each{ |v|
      replace = v[0].gsub(".","_")
      new_h[replace] = v[1]
      h.delete(v[0])
    }
    new_h.each{ |v|
      h[v[0]] = v[1]
    }
    @statements.each{ |s|
      if isStatement(s)
        lexp = get_name(s)
        type = nil
        h.each{ |v|
          iotype = get_iotype_from_hash(v)
          fdpsname = get_fdpsname_from_hash(v)
          if iotype == "FORCE" && fdpsname == lexp
            type = get_type_from_hash(v)
          end
        }
        s.replace_name(lexp,"FORCE_"+lexp) if type != nil

        append_iotag(s,"EPI",h)
        append_iotag(s,"EPJ",h)
        append_iotag(s,"FORCE",h)

      elsif s.class == TableDecl
        s.add_to_varhash
      end
    }
    #@statements.each{ |s|
    #  p s.convert_to_code("reference")
    #}
  end

  def process_iodecl_alias(ios,h=$varhash)
    ios.each{|x|
      if x.modifier == nil && ["EPI","EPJ","FORCE"].index(x.iotype)
        isIncluded = false
        h.each{ |v|
          fdpsname = v[1][2]
          isIncluded = true if fdpsname == x.fdpsname
        }
        abort "alias to undefined member #{x.fdpsname}" if !isIncluded
        h[x.name] = [x.iotype, x.type, x.fdpsname, "alias"]
      else
        h[x.name] = [x.iotype, x.type, x.fdpsname, x.modifier]
      end
    }
  end

end


