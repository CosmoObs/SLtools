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
	sed -i "s/%PROJECT_NAME%/$application/" doc/doxygen.cfg
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

[ ! -z "$1" ] && opt1="$1" || opt1="null"
[ ! -z "$2" ] && opt2="$2" || opt2="--all"

[ "$opt2" == "--all" -o "$opt2" == "--no-doc" -o "$opt2" == "--no-pdf" ] && mode=$2 || opt1="null"

if [ "$opt1" == "sltools" ]
then

    # Read package VERSION and exec building function:
    version=$(cat ./etc/version.txt)
    folder=$opt1-v$version

    echo "Building $opt1 v$version..."
    $1 $folder

    find $folder -name "*.pyc" -delete

    if [ "$mode" != "--no-doc" ]; then
	# compile documentation
	echo "Generating Doxygen documentation..."
	GenDoxygenDoc $opt1 $version $mode
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

    echo "Usage: ./build.sh { sltools } [ --all | --no-doc | --no-pdf ]"
    echo
    exit 1

fi


