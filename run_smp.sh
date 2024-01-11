#!/bin/bash
flag=$1

if [ ! $# == 1 ]; then
    echo "Usage: $0 OPTION={STANDALONE,WITHTOPMODEL}"
    echo "One of these options must be specified to run SMP"
  exit
fi

if [ $flag == "STANDALONE" ] || [ "$flag" == "WITHTOPMODEL" ]; then
echo "SMP running with option $flag"
else
echo "Invalid option! $flag"
exit
fi


args=" "
exe_name=" "
if [ $flag == "STANDALONE" ]; then
    args='./configs/config_conceptual.txt'
    exe_name='smp_standalone'
else if [ $flag == "WITHTOPMODEL" ]; then
	 args="configs/topmod.run configs/config_topmodel.txt"
	 exe_name='smp_topmodel'
     fi
fi
echo "config file: $args"
./build/${exe_name} $args
