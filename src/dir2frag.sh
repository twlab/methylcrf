#!/bin/bash
#------------------------------------------------------------------------#
# Copyright 2012                                                         #
# Author: stevens _at_ cse.wustl.edu (from Ting Wang original)           #
#                                                                        #
# This program is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation, either version 3 of the License, or      #
# (at your option) any later version.                                    #
#                                                                        #
# This program is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# GNU General Public License for more details.                           #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.  #
#------------------------------------------------------------------------#

# req 
#  olapBed

tmp=$(mktemp $0.XXXXXX.tmp)
cat $1/* |cut -f1,2,3|sort -u >$tmp
olapBed -c $tmp $2 |awk '{OFS="\t";print $4,!!$7}' 
rm $tmp
