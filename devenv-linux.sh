
GMSHSDK=~/local/gmsh-4.1.3-Linux64-sdk

export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:${PATH}/bin
export INCLUDE=${GMSHSDK}/include:${INCLUDE}
export LIB=${GMSHSDK}/lib:${LIB}
export PYTHONPATH=${GMSHSDK}/lib:${PYTHONPATH}

# examples in ${GMSHSDK}/share/doc/gmsh/demos/api/