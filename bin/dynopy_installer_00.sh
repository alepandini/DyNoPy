#!/bin/bash

#
# main dir for installation
#
dir_main="$(pwd)"
#
# hh-suite installation variables
#
dir_hh="hh-suite"
dir_hhbin="$dir_hh/bin"
dir_hhdbs="$dir_hh/database"
# change this variable to download a different database
file_dbs="UniRef30_2020_02_hhsuite"
#
# ccmpred variables
#
dir_ccm="$dir_main/CCMpred"
#
# cmake variables
#
dir_cmk="$dir_main/cmake-3.17.0"
echo "********************************"
echo "*                              *"
echo "* INSTALLING CMAKE             *"
echo "*            CCMpred           *"
echo "*            HH-BLITS          *"
echo "*                              *"
echo "********************************"
echo "USER WILL REQUIRE sudo ACCESS..."
echo "... TO INSTALL cmake            "
#
# install hh-suite based on instructions from https://github.com/soedinglab/hh-suite 
#
hhsuite_install()
{
	if [ ! -d $dir_hh ];
	then
		git clone https://github.com/soedinglab/hh-suite.git
		mkdir -p hh-suite/build && cd hh-suite/build
		cmake -DCMAKE_INSTALL_PREFIX=. ..
		make -j 4 && make install
		export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
	fi
	#
	# download the uniprot30 database into $dir_hhdbs
	# check if folder has content gt 100GB (ideal 12)
	#
	if [ -d $dir_hhdbs ]
	then
		cd $dir_hhdbs 
		#if [ ! -f "UniRef30_2020_02_cs219.ffdata" ];
		((folder_size=$(du -s | awk '{print $1}')))
		if [ $folder_size -gt 140034292 ];
		then
			echo "$file_dbs DATABASE DOWNLOAD INCOMPLETE. WILL TRY AGAIN"
			cd $dir_main
			rm -rf $dir_hhdbs
		fi
	fi
	if [ ! -d $dir_hhdbs ];
	then
		# create and enter the db dir
		mkdir $dir_hhdbs && cd $dir_hhdbs
		#delete the line below after debugging
		#rm $file_dbs
		if [ ! -f $file_dbs ]
		then
			wget http://wwwuser.gwdg.de/~compbiol/uniclust/2020_02/$file_dbs.tar.gz
			tar -xvf $file_dbs.tar.gz
			mkdir $file_dbs
			mv *.* $filde_dbs
		fi
		
	fi
}

ccm_install()
{
	if [ ! -d $dir_ccm ]
	then	
		git clone --recursive https://github.com/soedinglab/CCMpred.git
		cd $dir_ccm
		mkdir build && cd build
		cmake ../
		make
	fi
}
#
# installing ccmpred requires cmake 3.17
# need sudo permissions to run this phase
#
cmake_install()
{
	#rm -rf $dir_cmk
	# dependencies identified on meriadoc
	sudo make apt install libssl-dev
	file_cmk="cmake-3.17.0.tar.gz"
	if [ ! -f $file_cmk ];
	then
		wget https://github.com/Kitware/CMake/releases/download/v3.17.0/cmake-3.17.0.tar.gz
		tar -xvf cmake-3.17.0.tar.gz
	fi
	if [ -d $dir_cmk ];
	then

		cd $dir_cmk
		./configure
		sudo make install
	fi

}
#
# install order
# cmake (comment this function if you have cmake 3.17 need sudo access)
# ccmpred (dependent on cmake 3.17)
# hh-suite
#
cmake_install
cd $dir_main
ccm_install
cd $dir_main
hhsuite_install
