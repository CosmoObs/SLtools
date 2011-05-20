#!/bin/bash

############################
# Build script for sltools
#
# Angelo Fausti
# Carlos Brandt
############################

namespace2module()
{   # After have compiled documentation with doxygen,
    # modify some strings inside gtml/tex files.
    # File extension (html,tex) is passed through arg $1.

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

GenDoxygenDoc()
{   # Generate Doxygen documentation.
    # Receives app name in $1 (e.g., sltools), version number in $2,
    # and $3 gives the documentation 'mode' (pdf/web).

	project_name="SLtools"
	application=$1
	version=$2
	mode=$3

	# Check if doxypy is present 
	[ $(which doxypy.py) ] || { echo "Failed: doxypy not found."; exit 1; }

	dir=$(pwd)
	cd $FOLDER &> /dev/null
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
        cd doc/html && namespace2module html
	)

	(
        cd doc/latex && namespace2module tex

	    if [ "$mode" != "--no-pdf" ]
        then
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

	# Clone non-directory entries, but now make links instead of copying the files
	for i in `find ./[_a-zA-Z0-9]* -name "*"`
	do
		[ ! -d $i ] && ln -sf $PWD/$i $DEST/$i
	done

}

clone_devel()
{   # Clone $opt (e.g., sltools) tree to $destdir directory.

    opt=$1
    destdir=$2
    
    [ ! -w $destdir -o "$destdir" == "$PWD" ] && { echo "Directory $destdir is not writeable. Try again."; exit 1; }
    
    destdir="${destdir}/${opt}"
    mkdir -p $destdir

    Link_subtree_Clone $destdir

    rm -f $destdir/build.sh
    rm -f $destdir/*_TO_INCLUDE_IN_BUILD

}

build_sltools()
{   # SLtools building function
    # It receives the building $mode and the destination $FOLDER folder

    mode=$1
	FOLDER=$2

	if [ -d "$FOLDER" -o -f "${FOLDER}.tgz" ]
	then

		echo "Removing old build ..."
		rm -rf $FOLDER{,.t*gz} &> /dev/null

	fi

    echo "Creating package folder $FOLDER ..."
	mkdir $FOLDER

    if [ -f "C_MODULES_TO_INCLUDE_IN_BUILD" ]
    then

        echo "Cleaning C modules directories ..."
        for i in `grep -v "^#" C_MODULES_TO_INCLUDE_IN_BUILD`
        do
            cd $i && { make clean; cd - &> /dev/null; }
        done

    fi
    
	echo "Copying files ..."
    if [ -f "DIRECTORIES_TO_INCLUDE_IN_BUILD" ]
    then

        cp __init__.py $FOLDER
        cp -r `grep -v "^#" DIRECTORIES_TO_INCLUDE_IN_BUILD` $FOLDER
        cp -r -p bin $FOLDER
        cp etc/* ${FOLDER}/. &> /dev/null

    else

        cp -r -p * ${FOLDER}/.

    fi

    find $FOLDER -name "*.pyc" -delete

    if [ "$mode" != "--no-doc" ]
    then

		# compile documentation
		echo "Generating Doxygen documentation ..."
		GenDoxygenDoc "$opt" "$version" "$mode"

    fi

    version=${FOLDER#*'-'}
    sed "s/%VERSION%/$version/" $FOLDER/install.sh > $FOLDER/install.tmp
    mv $FOLDER/install.tmp $FOLDER/install.sh
    chmod +x $FOLDER/install.sh

}


####################################
##### Start build's main block #####
####################################

# Call a specific build function:
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
				echo "Bleh: wrong argument for --devel option. Try an existing directory ;)"
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
			echo "Give me an option to work with."
			opt="null"
		;;		
	esac
	shift
done
# If no argument is given (./build sltools), $mode will not be set.
# No problem if we are properly using the variable (within "$quotes")

if [ "$opt" != "sltools" ]
then
    echo "Usage:  ./build.sh { sltools } [ --devel <DESTDIR> | --no-doc | --no-pdf  | -o <FILENAME>]"
    echo ""
	echo " (default: all)     If no argument is given the building system run defaults: build entire package"
	echo " --devel <path>     clone subtree structure to DESTINY(path). Useful for developing outside Git"
	echo " --no-doc           No Doxygen documentation is compiled"
	echo " --no-pdf           No Doxygen PDF pages are generated, just HTML doc"
	echo " -o <filename>	  Name of output file"
	echo ""
   
    exit 2
fi

if [ "$mode" == "--devel" ]
then

    echo "Cloning devel tree into $destdir ..."
    clone_devel "$opt" "$destdir"
    echo "Done."    
    exit 0
fi

version=$(head -n1 ./etc/version.txt)
FOLDER="${opt}-v${version}"

echo "Building $opt version: $version ..."
build_sltools "$mode" "$FOLDER"

echo "Compressing $FOLDER ..."
tar czf "${FOLDER}.tgz" "$FOLDER"
rm -rf $FOLDER

echo "Done."
