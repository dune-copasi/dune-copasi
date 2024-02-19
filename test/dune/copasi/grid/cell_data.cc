
#include <dune/copasi/grid/cell_data.hh>

#include <dune/common/parallel/mpihelper.hh>

#include "fixture.hh"

template<class GV>
void test_cell_data(GV grid_view, bool do_shring_to_fit) {
  Dune::Copasi::CellData<GV, int> cell_data(grid_view);

  auto shrink_to_fit = [&]{
    if (do_shring_to_fit and cell_data.capacity() != cell_data.size()) {
      cell_data.shrink_to_fit();
      EXPECT_EQ(cell_data.size(), cell_data.capacity());
    }
  };

  std::vector<int> data(grid_view.size(0), 0);
  cell_data.addData("null", data);
  shrink_to_fit();

  std::iota(std::begin(data), std::end(data), 0);
  cell_data.addData("iota", data);
  shrink_to_fit();

  auto is_pair = [](int v){return (v % 2) == 0;};
  auto is_odd = [](int v){return (v % 2) != 0;};
  cell_data.addData("iota_pair", data, is_pair);
  shrink_to_fit();
  
  cell_data.addData("iota_odd", data, is_odd);
  shrink_to_fit();

  EXPECT_EQ(cell_data.size(), 4);

  std::vector<int> data_buffer;
  std::vector<bool> mask_buffer;
  std::size_t iota = 0;
  for (const auto& cell : elements(grid_view)) {

    EXPECT_EQ(cell_data.getData(0, cell).value(), 0);
    EXPECT_EQ(cell_data.getData("null", cell).value(), 0);
    EXPECT_EQ(cell_data.getData("iota", cell).value(), iota);
    EXPECT_EQ(cell_data.getData("iota_pair", cell).has_value(), is_pair(iota));
    EXPECT_EQ(cell_data.getData("iota_odd", cell).has_value(), is_odd(iota));
    if (is_pair(iota))
      EXPECT_EQ(cell_data.getData("iota_pair", cell).value(), iota);
    else
      EXPECT_EQ(cell_data.getData("iota_odd", cell).value(), iota);

    cell_data.getData(cell, data_buffer, mask_buffer);
    
    EXPECT_TRUE(mask_buffer[0]);
    EXPECT_TRUE(mask_buffer[1]);
    EXPECT_EQ(mask_buffer[2], is_pair(iota));
    EXPECT_EQ(mask_buffer[3], is_odd(iota));

    EXPECT_EQ(data_buffer[0], 0);
    EXPECT_EQ(data_buffer[1], iota);
    if (is_pair(iota))
      EXPECT_EQ(data_buffer[2], iota);
    else
      EXPECT_EQ(data_buffer[3], iota);

    ++iota;
  }

  // cell data of all kind of grid views must be able to respond to the leaf view
  for (const auto& cell : elements(grid_view.grid().leafGridView())) {
    EXPECT_EQ(cell_data.getData(0, cell).value(), 0);
    EXPECT_EQ(cell_data.getData("null", cell).value(), 0);

    iota = cell_data.getData("iota", cell).value();
    EXPECT_LT(iota, grid_view.size(0));
    EXPECT_EQ(cell_data.getData("iota", cell).value(), iota);
    EXPECT_EQ(cell_data.getData("iota_pair", cell).has_value(), is_pair(iota));
    EXPECT_EQ(cell_data.getData("iota_odd", cell).has_value(), is_odd(iota));
    if (is_pair(iota))
      EXPECT_EQ(cell_data.getData("iota_pair", cell).value(), iota);
    else
      EXPECT_EQ(cell_data.getData("iota_odd", cell).value(), iota);

    cell_data.getData(cell, data_buffer, mask_buffer);

    EXPECT_TRUE(mask_buffer[0]);
    EXPECT_TRUE(mask_buffer[1]);
    EXPECT_EQ(mask_buffer[2], is_pair(iota));
    EXPECT_EQ(mask_buffer[3], is_odd(iota));

    EXPECT_EQ(data_buffer[0], 0);
    EXPECT_EQ(data_buffer[1], iota);
    if (is_pair(iota))
      EXPECT_EQ(data_buffer[2], iota);
    else
      EXPECT_EQ(data_buffer[3], iota);
  }
}

// declare a dummy class that forwards it's template argument as the fixture to test
template <typename Fixture> class ForwardFixture : public Fixture {};

// declare the templated fixture
TYPED_TEST_SUITE_P(ForwardFixture);

// declare test to make on each space
TYPED_TEST_P(ForwardFixture, TestCellData) {
  for (auto do_shring_to_fit : {false, true}) {
    test_cell_data(this->_grid->leafGridView(), do_shring_to_fit);
    for (std::size_t l = 0; l != this->_grid->maxLevel(); ++l)
      test_cell_data(this->_grid->levelGridView(l), do_shring_to_fit);
  }
}


// register the tests
REGISTER_TYPED_TEST_SUITE_P(ForwardFixture, TestCellData);

// instantiate the test for each type
INSTANTIATE_TYPED_TEST_SUITE_P(CellData, ForwardFixture, GridFixtures);

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  Dune::MPIHelper::instance(argc, argv);
  return RUN_ALL_TESTS();
}
