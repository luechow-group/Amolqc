" Vim syntax file
" Language: Amolqc input file
" Maintainer: Leonard Reuter
" Latest Revision: 14 May 2020


if exists("b:current_syntax")
  finish
endif

syn match       amiDefault "\S"

syn keyword	amiTodo	FIXME NOTE NOTES TODO XXX YYY ZZZ contained
syn match	amiComment	"!.*$" contains=amiTodo

" Logicals. Note, that fortran also accepts other strings (eg. 'f', 't')
syn match amiLogical '\.True\.\|\.true\.\|\.TRUE\.\|\.T\.\|\.t\.' contained
syn match amiLogical '\.False\.\|\.false\.\|\.FALSE\.\|\.F\.\|\.f\.' contained

" Integer with -,+ or nothing in front and k, M or nothing in back
syn match amiNumber '[+-]\{0,1}\d\+[kM]\{0,1}' contained

" Floating point number with decimal, without e,E,d,D 
syn match amiNumber '[+-]\{0,1}\d\+\.\d*' contained
syn match amiNumber '[+-]\{0,1}\d*\.\d\+' contained

" Floating point number without decimal, with and e,E,d,D
syn match amiNumber '[+-]\{0,1}\d\+[EeDd][+-]\{0,1}\d\+' contained

" Floating point number with decimal, with e,E,d,D
syn match amiNumber '[+-]\{0,1}\d\+\.\d*[EeDd][+-]\{0,1}\d\+' contained
syn match amiNumber '[+-]\{0,1}\d*\.\d\+[EeDd][+-]\{0,1}\d\+' contained

" Strings
syn region amiString start='"' end='"' contains=amiTodo contained
syn region amiString start="'" end="'" contains=amiTodo contained

" Function Calls
syn region amiFunction start='\$' end='('me=s-1

" Macro Calls, includes selfmade macros startin with 'macro'
syn region amiMacro start='\$generate_sample' end='('me=s-1
syn region amiMacro start='\$generate_walker' end='('me=s-1
syn region amiMacro start='\$jastrow_emin_lin' end='('me=s-1
syn region amiMacro start='\$jastrow_emin_lm' end='('me=s-1
syn region amiMacro start='\$jastrow_emin_nw' end='('me=s-1
syn region amiMacro start='\$jastrow_varmin_fast' end='('me=s-1
syn region amiMacro start='\$jastrow_varmin_safe' end='('me=s-1
syn region amiMacro start='\$jas_emin_lin' end='('me=s-1
syn region amiMacro start='\$jas_emin_lm' end='('me=s-1
syn region amiMacro start='\$jas_emin_nr' end='('me=s-1
syn region amiMacro start='\$jas_emin_snr' end='('me=s-1
syn region amiMacro start='\$jas_varmin_fast' end='('me=s-1
syn region amiMacro start='\$jas_varmin_safe' end='('me=s-1
syn region amiMacro start='\$macro' end='('me=s-1

" Key and Value within function calls
syn region amiKey start='[(,]'ms=e+1 end='[=,)]'me=s-1 contains=amiComment contained
syn region amiValue start='='ms=e+1 end='[,)]'me=s-1 contains=amiString,amiNumber,amiLogical,amiComment contained
syn region amiInFunction start='(' end=')' contains=amiKey,amiValue transparent

let b:current_syntax = "amolqcami"

hi def link amiTodo	Todo
hi def link amiComment	Comment
hi def link amiNumber	Constant
hi def link amiString	Special
hi def link amiFunction	Function
hi def link amiKey	Keyword
hi def link amiValue    Type
hi def link amiDefault	Error
hi def link amiLogical	PreProc
hi def link amiMacro	Underlined

