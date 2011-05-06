#!/bin/bash

# Build script for sltools

# Post process Doxygen documentation,
# replace "namespace/package" ocurrences with "module"
#
# Used by GenDoxygenDoc()

namespace2module()
{

	ext=$1

	sed -i 's/Namespace/Module/' *.$ext &> /dev/null 
	sed -i 's/Package/Module/' *.$ext &> /dev/null
	sed -i 's/packages/modules/' *.$ext &> /dev/null
	sed -i 's/namespace members/module members/' *.$ext &> /dev/null
	sed -i 's/namespaces they/modules they/' *.$ext &> /dev/null

	if [ "$1" = "html" ]; then
		sed -i 's/package/module/' *.$ext &> /dev/null
	fi

	if [ "$1" = "tex" ]; then
		sed -i 's/^package/module/' *.$ext &> /dev/null
	fi
}

# Generate Doxygen documentation

GenDoxygenDoc()
{
	project_name="SLtools"
	application=$1
	version=$2
	mode=$3

	# Check if doxypy is present 
	[ $(which doxypy.py) ] || { echo "Failed: doxypy not found."; exit 1; }

	dir=$(pwd)
	cd $folder &> /dev/null
	mkdir doc
	cp ../doc/doxygen.cfg doc/

	# Edit doxygen configuration file for the current application
	sed -i "s/%PROJECT_NAME%/$project_name/" doc/doxygen.cfg
	sed -i "s/%PROJECT_NUMBER%/\"Version $version\"/" doc/doxygen.cfg


	# work around for bug #600 [AddArcs] Fix doxygen documentation issue

	for file in $(find . -name "__init__.py")
	do
		mv $file $(dirname $file)/__init__
	done

	
	doxygen doc/doxygen.cfg &> /dev/null || { echo "Failed.";  cd $dir &> /dev/null; exit 1; }

	find . -name "__init__" -exec mv '{}' '{}'.py \;

	(
	    cd doc/html
	    namespace2module html
	)

	(
        cd doc/latex &> /dev/null
	    namespace2module tex

	    if [ "$mode" != "--no-pdf" ]; then
		make &> /dev/null || { echo "Failed to create PDF documentation."; exit 1; }
		[ -f refman.pdf ] && mv refman.pdf $application-v$version.pdf
	    fi
	)

	cd $dir &> /dev/null

}

Link_subtree_Clone()
{
	DEST=$1
	
	# Clone directory entries below current dir to dir passed trough argument '$1'
	for i in `find ./[a-zA-Z0-9]* -name "*"`
	do
		[ -d $i ] && mkdir -p $DEST/$i
	done

	# Clone non-directory entries, but now make link instead of copying the files
	for i in `find ./[_a-zA-Z0-9]* -name "*"`
	do
		[ ! -d $i ] && ln -sf $PWD/$i $DEST/$i
	done

}

# SLtools build

sltools()
{

	folder=$1

	if [ -d "$folder" -o -f "$folder.tgz" ]
	then
		echo "Removing old build..."
		rm -rf $folder* &> /dev/null
	fi

	mkdir $folder

	echo "Creating missing sltools folders..."

	mkdir -p $folder

        cd lens/get_nfw_concentration &> /dev/null
	make clean &> /dev/null
	cd - &> /dev/null

	cd lens/compute_nfw_lens_parameters &> /dev/null
	make clean &> /dev/null
	cd - &> /dev/null

	echo "Copying files..."

	cp __init__.py $folder
	cp -r `grep -v "^#" DIRECTORIES_TO_INCLUDE_IN_BUILD` $folder  
	cp -r -p bin $folder
	cp -r etc/install.sh $folder
	cp -r etc/INSTALL $folder
	cp -r etc/README $folder
	cp -r etc/HISTORY $folder

}

# Call a specific build function

while test -n "$1"
do
	case "$1" in
		--no-doc)
			mode="$1"
		;;
		--no-pdf)
			mode="$1"
		;;
		--all)
			mode="$1"
		;;
		sltools)
			opt="sltools"
		;;
		--devel)
			mode="$1"
			shift
			destdir="$1"
			if [ ! -d "$destdir" -o -z "$destdir" ]
			then
				echo "Bleh: wrong argument for --devel option."
				exit 1
			fi
		;;
		-o)
			shift
			folder="$1"
		;;
		-h | --help)
			opt="null"
		;;
		*)
			echo "Give an option to work with."
			opt="null"
		;;		
	esac
	shift
done

###
#[ ! -z "$1" ] && opt1="$1" || opt1="null"
#[ ! -z "$2" ] && opt2="$2" || opt2="--all"
#[ ! -z "$3" ] && opt3="$3" || opt3=""

#[ "$opt2" == "--devel" -o "$opt2" == "--no-doc" -o "$opt2" == "--no-pdf" -o "$opt2" == "--all" ] && mode=$2 || opt1="null"
#[ "$opt2" == "--devel" -a -z "$opt3" ] && opt1="null" || destdir=$3

if [ "$opt" == "sltools" ]
then

	if [ "$mode" == "--devel" ]; then
		[ ! -w $destdir -o "$destdir" == "$PWD" ] && { echo "Directory $destdir is not writeable. Try again."; exit 1; }
		echo "Cloning devel tree into $destdir ..."
		mkdir -p $destdir/sltools
		destdir="$destdir/sltools"
		Link_subtree_Clone $destdir
		rm -f $destdir/build.sh
		rm -f $destdir/DIRECTORIES_TO_INCLUDE_IN_BUILD
		echo "Done."
		exit 0
	fi
	
    # Read package VERSION and exec building function:
    version=$(cat ./etc/version.txt)
    [ -z "$folder" ] && folder=$opt-v$version

    echo "Building $opt version: $version..."
    sltools $folder

    find $folder -name "*.pyc" -delete

    if [ "$mode" != "--no-doc" ]; then
		# compile documentation
		echo "Generating Doxygen documentation..."
		GenDoxygenDoc $opt $version $mode
    fi

    sed "s/%VERSION%/$version/" $folder/install.sh > $folder/install.tmp
    mv $folder/install.tmp $folder/install.sh
    chmod +x $folder/install.sh

    echo "Compressing $folder..."
    tar cvzf $folder.tgz $folder &> /dev/null
    rm -rf $folder &> /dev/null

    echo "Done."
    echo

else

    echo "Usage: ./build.sh { sltools } [ --devel <DESTDIR> | --no-doc | --no-pdf  | -o <FILENAME>]"
    echo ""
	echo " (default: all)     If no argument is given the building system run defaults: build entire package"
	echo " --devel <path>     clone subtree structure to DESTINY(path). Useful for developing outside Git"
	echo " --no-doc           No Doxygen documentation is compiled"
	echo " --no-pdf           No Doxygen PDF pages are generated, just HTML doc"
	echo " -o <filename>	  Name of output file"
	echo ""
    exit 1

fi


