# Includes the Package, bringing the lx_functions into scope
using FranklinTheorems

# Includes the custom markdown files, bringing the `\newcommand` and `\newenvironment` definitions into scope.
Franklin.include_external_config(FranklinTheorems.config_path())

# LaTeX Preamble
Franklin.convert_md(read("./_fkl_libs/preamble.texmd", String); isconfig=true)

function hfun_bar(vname)
  val = Meta.parse(vname[1])
  return round(sqrt(val), digits=2)
end

function hfun_m1fill(vname)
  var = vname[1]
  return pagevar("index", var)
end

function lx_baz(com, _)
  # keep this first line
  brace_content = Franklin.content(com.braces[1]) # input string
  # do whatever you want here
  return uppercase(brace_content)
end
