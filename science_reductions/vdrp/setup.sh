# allow python to find vdrp modules
export PYTHONPATH=$PYTHONPATH:/work/00115/gebhardt/maverick/vdrp/vdrp/src/python
# this is required for daophot
#module load gcc/5.4.0
# to find daophot
export PATH=~gebhardt/lib/daoprogs:$PATH
# to find daomaster
export PATH=~gebhardt/lib/daoprogs/moreprogs2:$PATH
# to find getoff2
export PATH=/work/00115/gebhardt/maverick/sci/progs:$PATH
# to find  immosaicv & imrot
export PATH=~gebhardt/bin:$PATH
