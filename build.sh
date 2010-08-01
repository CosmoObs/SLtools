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

	# compile documentation
	echo "Generating Doxygen documentation for $folder"

	# Check if doxypy is present 
	[ $(which doxypy.py) ] || { echo "Failed: doxypy not found."; exit 1; }

	dir=$(pwd)
	cd $folder &> /dev/null
	mkdir doc
	cp ../doc/doxygen.cfg doc/

	# Edit doxygen configuration file for the current application
	sed -i "s/%PROJECT_NAME%/$application/" doc/doxygen.cfg
	sed -i "s/%PROJECT_NUMBER%/\"Version $version\"/" doc/doxygen.cfg


	GenDoxygenDoc sltools $version


	# work around for bug #600 [AddArcs] Fix doxygen documentation issue

	for file in $(find . -name "__init__.py")
	do
		mv $file $(dirname $file)/__init__
	done

	
	doxygen doc/doxygen.cfg &> /dev/null || { echo "Failed.";  cd $dir &> /dev/null; exit 1; }

	find . -name "__init__" -exec mv '{}' '{}'.py \;	

	cd doc/html
	namespace2module html

	cd ../latex &> /dev/null
	namespace2module tex 

	make &> /dev/null || { echo "Failed to create PDF documentation."; exit 1; }
	[ -f refman.pdf ] && mv refman.pdf $application-v$version.pdf
	cd $dir &> /dev/null

}

# SLtools build

sltools()
{

	version=$1

	echo "Building sltools v$version"

	folder=sltools-v$version

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
        cp -r catalog $folder  
        cp -r gravlens $folder  
        cp -r io $folder	
        cp -r lens $folder
        cp -r plot $folder
        cp -r string $folder
        cp -r image $folder
        cp -r coordinate $folder
        cp -r bin $folder
	cp -r etc/install.sh $folder
	cp -r etc/INSTALL $folder
	cp -r etc/README $folder
	cp -r etc/HISTORY $folder


	# Edit installation file version number

	sed -i "s/%VERSION%/$version/" $folder/install.sh

	# make package
	echo "Compressing $folder..."
	tar cvzf $folder.tgz $folder &> /dev/null
	rm -rf $folder &> /dev/null

	echo "Done."
	echo
}


# Call a specific build function

if [ "$1" == "sltools" ]
then

    # Read package VERSION and exec building function:

    version=$(cat ./etc/version.txt)

    $1 $version

else

    echo "Usage: ./build.sh [ sltools ]"
    echo
    exit 1

fi


