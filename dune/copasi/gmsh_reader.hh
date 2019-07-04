#ifndef DUNE_COPASI_GMSH_READER_HH
#define DUNE_COPASI_GMSH_READER_HH

#include <dune/grid/multidomaingrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/common/gridinfo.hh>

#include <type_traits>
#include <algorithm>

namespace Dune::Copasi {

template<class Grid>
class GmshReader;

template<class HostGrid, class MDGTraits>
class GmshReader<Dune::mdgrid::MultiDomainGrid<HostGrid,MDGTraits>>
{
public:
  using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid,MDGTraits>;

  static std::pair<Grid*,HostGrid*> read (const std::string& fileName, bool verbose = true, bool insertBoundarySegments=true)
  {
    // make a grid factory
    Dune::GridFactory<HostGrid> factory;

    // create parse object
    GmshReaderParser<HostGrid> parser(factory,verbose,insertBoundarySegments);
    parser.read(fileName);

    auto index_map = parser.elementIndexMap();

    HostGrid* host_grid = factory.createGrid();

    int max_subdomains = 1+*std::max_element(index_map.begin(),index_map.end());
    MDGTraits* traits;
    if constexpr (std::is_default_constructible_v<MDGTraits>)
      traits = new MDGTraits;
    else
      traits = new MDGTraits(max_subdomains);

    Grid* grid = new Grid(*host_grid,*traits);
    delete traits;

    grid->startSubDomainMarking();
    uint i = 0;
    for (const auto& cell : elements(grid->leafGridView()))
      grid->addToSubDomain(index_map[i]-1,cell), i++;
    
    grid->preUpdateSubDomains();
    grid->updateSubDomains();
    grid->postUpdateSubDomains();

    for (int i = 0; i < max_subdomains; ++i)
      gridinfo(grid->subDomain(i));

    return std::make_pair(grid,host_grid);
  }















  // /** \todo doc me */
  // static Grid* read (const std::string& fileName,
  //                    std::vector<int>& boundarySegmentToPhysicalEntity,
  //                    std::vector<int>& elementToPhysicalEntity,
  //                    bool verbose = true, bool insertBoundarySegments=true)
  // {
  //   // make a grid factory
  //   Dune::GridFactory<Grid> factory;

  //   // create parse object
  //   GmshReaderParser<Grid> parser(factory,verbose,insertBoundarySegments);
  //   parser.read(fileName);



  //   boundarySegmentToPhysicalEntity.swap(parser.boundaryIdMap());
  //   elementToPhysicalEntity.swap(parser.elementIndexMap());

  //   return factory.createGrid();
  // }

  // /** \todo doc me */
  // static void read (Dune::GridFactory<Grid>& factory, const std::string& fileName,
  //                   bool verbose = true, bool insertBoundarySegments=true)
  // {
  //   // create parse object
  //   GmshReaderParser<Grid> parser(factory,verbose,insertBoundarySegments);
  //   parser.read(fileName);
  // }

  // /** \todo doc me */
  // static void read (Dune::GridFactory<Grid>& factory,
  //                   const std::string& fileName,
  //                   std::vector<int>& boundarySegmentToPhysicalEntity,
  //                   std::vector<int>& elementToPhysicalEntity,
  //                   bool verbose = true, bool insertBoundarySegments=true)
  // {
  //   // create parse object
  //   GmshReaderParser<Grid> parser(factory,verbose,insertBoundarySegments);
  //   parser.read(fileName);

  //   boundarySegmentToPhysicalEntity.swap(parser.boundaryIdMap());
  //   elementToPhysicalEntity.swap(parser.elementIndexMap());
  // }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GMSH_READER_HH