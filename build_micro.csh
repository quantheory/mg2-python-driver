#!/usr/bin/env csh

#This script uses f2py to compile the MG2 microphysics code and 
#sets up the necessary linkage to be able to call it from python.
#
#Problems: 
#1. ACME modules tend to define "r8 = selected_real_kind(12)" near their top, then 
#define all subsequent variables using "real(r8) <varname>". This causes 2 problems:
#
#a. f2py creates a temporary f90 file which just lists the interface of each subroutine/function 
#you want to wrap. This ends up including lots of real(r8) mentions because most of the input/output 
#variables are r8 type... but r8 isn't defined in this file because the file only includes input/output 
#variables. The easiest way to fix this is to create a .pyf file, then search/replace real(kind=r8) with 
#real(kind=8) everywhere in the file. This is nice because it avoids changing the actual source code 
#you're trying to evaluate. Note, however, that it could screw up the intended precision on some machines/compilers... so be careful!
#
#b. f2py doesn't like the selected_real_kind function because it needs to map 
#f90 to c data types and having this defined through a function is too complicated. So
#you need to define the typemap (see https://sysbio.ioc.ee/projects/f2py2e/FAQ.html#q-what-if-fortran-90-code-uses-type-spec-kind-kind) for discussion. I added the line dict(real=dict(r8="double")) to .f2py_f2cmap in
#the directory I'm calling f2py from.
#
#2. By default, f2py tries to wrap both public and private subroutines in all f90 files you give it. 
#Private subroutines cause errors because the temporary f90 file f2py creates tries to "use" these 
#private subroutines. Thus it is best to explicitly specify the subroutines you want to wrap and to stick
#to public ones. Otherwise, modify the code to make the private subroutines public. Note that you need to 
#use spaces rather than commas to specify a series of subroutines to wrap... commas make the code segfault 
#upon module load!


#BASIC SET UP:
#=================
# NM will be the name of the module. In python, subroutines will be 
# accessible under <NM>.<filename>.<fortran subroutine name>.

set NM=mg2

#This is the source code directory.
set DIR=.

#These are the files to load:
set FILES="${DIR}/wv_sat_methods.F90 ${DIR}/micro_mg_utils.F90 ${DIR}/micro_mg2_0.F90"

#NOTE: The line below defines all the functions/subroutines we want python 
#access to. Each of these functions must be defined as 'public' in the 
#fortran code we are parsing. Note also that the functions/subroutines below 
#are separated by spaces - the code still compiles if commas are used instead,
#but it crashes when you try to load it into python. 

set PUB_FNS='wv_sat_methods_init micro_mg_init micro_mg_tend calc_precip_frac micro_mg_utils_init'

#CREATE .pyf FILE WITH INTERFACES TO PYTHON:
#=================
#The .pyf file is a fortran file which just contains the interfaces needed
#to link the PUB_FNS to python.

f2py ${FILES} -m ${NM} --overwrite-signature -h ${NM}.pyf only: $PUB_FNS

#MODIFY .pyf FILE TO GET RID OF r8 REFERENCES:
#=================
#Because the .pyf file only contains info specific to interfaces, it omits the 
#"integer, parameter :: r8 = selected_real_kind(12)" definition at the beginning
#of most acme .F90 files. This causes compilation to fail because all the 
#variables in the interface are defined as "real(r8)". As a fix, I change
#all of these mentions in the .pyf file to "kind=8". This should work because 
#r8 is supposed to be double precision, which is what kind=8 means... but it 
#isn't as platform independent as selected_real_kind(12) so could screw up 
#python->fortran conversion on weird platforms. Still, this change is the only
#way I can think of reading the original (unmodified) code from ACME master into
#python and it should be obvious if type conversion is screwed up (because all 
#tests will fail) so I think it's ok.

sed -ie 's/kind=r8/kind=8/g' ${NM}.pyf

#On a similar note, sometimes the .pyf file truncates 'selected_real_kind(12)' to 
#'selected'. Replace these with kind=8 as well. 
sed -ie 's/real(kind=selected)/kind=8/g' ${NM}.pyf

#MODIFY .pyf FILE TO GET RID OF micro_mg_utils.F90 interface
#=================
#Some functions in micro_mg_utils.F90 (size_dist_param_liq, size_dist_param_basic, 
#and size_dist_param_ice) pass derived datatypes as input or output. f2py can't handle 
#derived datatypes and dies. The simple solution is to just omit these functions from 
#$PUB_FNS. 
#
#There are still problems because f2py tries to make all public variables in 
#a module available from python as well, which means it tries to make variables using
#the derived datatype available. This requires using sed to strip these references
#from the .pyf file. 

#HEAVY-HANDED SOLUTION: DELETE ENTIRE micro_mg_utils INTERFACE
#sed -i '/module micro_mg_utils/,/module micro_mg_utils/d' ${NM}.pyf

#JUST REMOVE LINES INCLUDING THE OFFENSIVE TYPE.
sed -i '/type(mghydrometeorprops)/d' ${NM}.pyf
sed -i '/type unknown_type/,/type unknown_type/d' ${NM}.pyf

#NOW ACTUALLY COMPILE THE CODE
#=================
#Note: -DHAVE_GAMMA_INTRINSICS #defines HAVE_GAMMA_INTRINSICS, which is 
#true for ifort and prevents an additional awkward dependency.

f2py -c --fcompiler=gnu95 -DHAVE_GAMMA_INTRINSICS ${NM}.pyf ${FILES}
