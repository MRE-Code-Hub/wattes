# WATTES - Wind And Tidal Turbine Embedded Simulator
# An actuator disc/line turbine model for coupling with CFD software
#
# Copyright (C) 2017 Heriot Watt University and the University of Edinburgh.
#
# Please see the AUTHORS file in the main source directory for a full list
# of copyright holders.
#
#	  Dr. A Creech
#	  Institute of Energy Systems
#	  School of Engineering
#	  University of Edinburgh
#	  
#	  angus_creech@hotmail.com
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
# -----------------------------------------------------------------------------


if [ "$1" == "" ]
then
    ncpus=4
else
    ncpus=$1
fi

rm -f output-*.log *~ output.log

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../lib


echo "run" > gdb.cmd

mpirun -np $ncpus xterm -sb -g 80x40 -e gdb -x gdb.cmd -se ../testprog > output.log

rm -f gdb.cmd

# mpirun -np $ncpus ../testprog
