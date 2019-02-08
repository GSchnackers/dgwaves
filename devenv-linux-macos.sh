
GMSHSDK=~/local/gmsh-sdk

export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:${PATH}
export INCLUDE=${GMSHSDK}/include:${INCLUDE}
export LIB=${GMSHSDK}/lib:${LIB}
export PYTHONPATH=${GMSHSDK}/lib:${PYTHONPATH}
export DYLD_LIBRARY_PATH=${GMSHSDK}/lib:${DYLD_LIBRARY_PATH}

# examples in ${GMSHSDK}/share/doc/gmsh/demos/api/
