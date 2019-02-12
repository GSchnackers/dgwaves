#include <cstdio>
#include <iostream>
#include <gmsh.h>

// [DOC from gmsh.h]
// A geometrical entity in the Gmsh API is represented by two integers: its
// dimension (dim = 0, 1, 2 or 3) and its tag (its unique, strictly positive
// identifier). When dealing with multiple geometrical entities of possibly
// different dimensions, the entities are packed as a vector of (dim, tag)
// integer pairs.
// typedef std::vector<std::pair<int, int> > vectorpair;

// [RB] note that entities of different dimensions can share the same tag
//      (although the pair (dim,tag) should be unique)

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 0;
    }

    // [DOC from gmsh.h]
    // Initialize Gmsh. This must be called before any call to the other functions in
    // the API. If `argc' and `argv' (or just `argv' in Python or Julia) are
    // provided, they will be handled in the same way as the command line arguments
    // in the Gmsh app. If `readConfigFiles' is set, read system Gmsh configuration
    // files (gmshrc and gmsh-options).
    // GMSH_API void initialize(int argc = 0, char ** argv = 0,
    //                          const bool readConfigFiles = true);
    gmsh::initialize(argc, argv);
    gmsh::option::setNumber("General.Terminal", 1); // enables "gmsh::logger::write(...)"
    gmsh::open(argv[1]);                            // reads the msh file

    // explore the mesh: what type of 2D elements do we have?

    // [DOC from gmsh.h]
    // Get the types of elements in the entity of dimension `dim' and tag `tag'.
    // If `tag' < 0, get the types for all entities of dimension `dim'. If `dim'
    // and `tag' are negative, get all the types in the mesh.
    // GMSH_API void getElementTypes(std::vector<int> & elementTypes,
    //                               const int dim = -1,
    //                               const int tag = -1);
    std::vector<int> eleTypes;
    gmsh::model::mesh::getElementTypes(eleTypes, 2);
    if (eleTypes.size() != 1)
    {
        gmsh::logger::write("Hybrid meshes not handled in this example!",
                            "error");
        // [RB] I guess that "getElementEdgeNodes" is not implemented for hybrid meshes.
        return 1;
    }

    // get all the properties (name, order, etc) related to this type of element

    // [DOC from gmsh.h]
    // Get the properties of an element of type `elementType': its name
    // (`elementName'), dimension (`dim'), order (`order'), number of nodes
    // (`numNodes') and parametric node coordinates (`parametricCoord' vector, of
    // length `dim' times `numNodes').
    // GMSH_API void getElementProperties(const int elementType,
    //                                    std::string & elementName,
    //                                    int & dim,
    //                                    int & order,
    //                                    int & numNodes,
    //                                    std::vector<double> & parametricCoord);
    int eleType2D = eleTypes[0];
    std::string name;
    int dim, order, numNodes;
    std::vector<double> paramCoord;
    gmsh::model::mesh::getElementProperties(eleType2D, name, dim, order,
                                            numNodes, paramCoord);
    gmsh::logger::write("2D elements are of type '" + name + "' (type = " +
                        std::to_string(eleType2D) + ") ");

    // iterate over all surfaces (i.e. entities of dim "2"), get the 2D elements
    // and create new 1D elements for all edges

    // [DOC from gmsh.h]
    // Get all the (elementary) geometrical entities in the current model. If `dim'
    // is >= 0, return only the entities of the specified dimension (e.g. points if
    // `dim' == 0). The entities are returned as a vector of (dim, tag) integer
    // pairs.
    // GMSH_API void getEntities(gmsh::vectorpair & dimTags,
    //                           const int dim = -1);
    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities, 2);
    std::cout << "There are " << entities.size() << " entitie(s) of dim 2\n";

    for (std::size_t i = 0; i < entities.size(); i++)
    {
        // [DOC from gmsh.h]
        // Get the elements of type `elementType' classified on the entity of of tag
        // `tag'. If `tag' < 0, get the elements for all entities. `elementTags' is a
        // vector containing the tags (unique, strictly positive identifiers) of the
        // elements of the corresponding type. `nodeTags' is a vector of length equal
        // to the number of elements of the given type times the number N of nodes
        // for this type of element, that contains the node tags of all the elements
        // of the given type, concatenated: [e1n1, e1n2, ..., e1nN, e2n1, ...]. If
        // `numTasks' > 1, only compute and return the part of the data indexed by
        // `task'.
        // GMSH_API void getElementsByType(const int elementType,
        //                                std::vector<int> & elementTags,
        //                                std::vector<int> & nodeTags,
        //                                const int tag = -1,
        //                                const size_t task = 0,
        //                                const size_t numTasks = 1);

        int s = entities[i].second;
        std::cout << "processing entity with tag #" << s << '\n';
        std::vector<int> elementTags, nodeTags;
        gmsh::model::mesh::getElementsByType(eleType2D, elementTags, nodeTags, s);
        gmsh::logger::write("- " + std::to_string(elementTags.size()) +
                            " elements in surface " + std::to_string(s));

        // get the nodes on the edges of the 2D elements

        // [DOC from gmsh.h]
        // Get the nodes on the edges of all elements of type `elementType'
        // classified on the entity of tag `tag'. If `primary' is set, only the
        // primary (begin/end) nodes of the edges are returned. If `tag' < 0, get the
        // edge nodes for all entities. If `numTasks' > 1, only compute and return
        // the part of the data indexed by `task'.
        // GMSH_API void getElementEdgeNodes(const int elementType,
        //                                   std::vector<int> & nodes,
        //                                   const int tag = -1,
        //                                   const bool primary = false,
        //                                   const size_t task = 0,
        //                                   const size_t numTasks = 1);

        std::vector<int> nodes;
        gmsh::model::mesh::getElementEdgeNodes(eleType2D, nodes, s);




        // create a new discrete entity of dimension 1

        // [DOC from gmsh.h]
        // Add a discrete geometrical entity (defined by a mesh) of dimension `dim' in
        // the current model. Return the tag of the new discrete entity, equal to `tag'
        // if `tag' is positive, or a new tag if `tag' < 0. `boundary' specifies the
        // tags of the entities on the boundary of the discrete entity, if any.
        // Specyfing `boundary' allows Gmsh to construct the topology of the overall
        // model.
        // GMSH_API int addDiscreteEntity(const int dim,
        //                                const int tag = -1,
        //                                const std::vector<int> & boundary = std::vector<int>());

        int c = gmsh::model::addDiscreteEntity(1);
        std::cout << "creating new discrete entity of dimension 1 (new tag=" << c << ")\n";

        // and add new 1D elements to it, for all edges

        // [DOC from gmsh.h]
        // Set the elements of type `elementType' in the entity of dimension `dim'
        // and tag `tag'. `elementTags' contains the tags (unique, strictly positive
        // identifiers) of the elements of the corresponding type. `nodeTags' is a
        // vector of length equal to the number of elements times the number N of
        // nodes per element, that contains the node tags of all the elements,
        // concatenated: [e1n1, e1n2, ..., e1nN, e2n1, ...]. If the `elementTag'
        // vector is empty, new tags are automatically assigned to the elements.
        // GMSH_API void setElementsByType(const int dim,
        //                                 const int tag,
        //                                 const int elementType,
        //                                 const std::vector<int> & elementTags,
        //                                 const std::vector<int> & nodeTags);

        int eleType1D = gmsh::model::mesh::getElementType("line", order);
        gmsh::model::mesh::setElementsByType(1, c, eleType1D, {}, nodes);

        // here we created two 1D elements for each edge; to create unique elements
        // it would be useful to call getElementEdgeNodes() with the extra `primary'
        // argument set to 'true' (to only get start/end nodes even in the
        // high-order case, i.e. consider topological edges), then sort them and
        // make them unique.

        // this could be enriched with additional info: each topological edge could
        // be associated with the tag of its parent element; in the sorting process
        // eliminating duplicates a second tag can be associated for internal edges,
        // allowing to keep track of neighbors
    }

    std::cout << "iterate over all 1D elements and get integration information\n";

    //gmsh::write("edges.msh");

    // iterate over all 1D elements and get integration information
    gmsh::model::mesh::getElementTypes(eleTypes, 1);

    // [DOC from gmsh.h]
    // Get the basis functions of the element of type `elementType' for the given
    // `integrationType' integration rule (e.g. "Gauss4") and `functionSpaceType'
    // function space (e.g. "IsoParametric"). `integrationPoints' contains the
    // parametric coordinates u, v, w and the weight q for each integeration
    // point, concatenated: [g1u, g1v, g1w, g1q, g2u, ...]. `numComponents'
    // returns the number C of components of a basis function. `basisFunctions'
    // contains the evaluation of the basis functions at the integration points:
    // [g1f1, ..., g1fC, g2f1, ...].
    // GMSH_API void getBasisFunctions(const int elementType,
    //                                 const std::string & integrationType,
    //                                 const std::string & functionSpaceType,
    //                                 std::vector<double> & integrationPoints,
    //                                 int & numComponents,
    //                                 std::vector<double> & basisFunctions);
    int eleType1D = eleTypes[0];
    std::vector<double> intpts, bf;
    int numComp;
    gmsh::model::mesh::getBasisFunctions(eleType1D, "Gauss3", "IsoParametric",
                                         intpts, numComp, bf);
    gmsh::model::getEntities(entities, 1);
    for (std::size_t i = 0; i < entities.size(); i++)
    {
        int c = entities[i].second;
        std::vector<int> elementTags, nodeTags;
        gmsh::model::mesh::getElementsByType(eleType1D, elementTags, nodeTags, c);
        gmsh::logger::write("- " + std::to_string(elementTags.size()) +
                            " elements on curve " + std::to_string(c));
        std::vector<double> jac, det, pts;
        gmsh::model::mesh::getJacobians(eleType1D, "Gauss3", jac, det, pts, c);
    }

    //gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}
