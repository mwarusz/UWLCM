#pragma once

#include <blitz/array.h>
#include <H5Cpp.h>
#include <map>

class PlotterCommon
{
  public:
  const string file;
  std::map<std::string, double> map;
  blitz::Array<float, 1> timesteps;
  double CellVol, DomainSurf;

  protected:
  H5::H5File h5f;
  H5::DataSet h5d;
  H5::Group h5g;
  H5::DataSpace h5s;

  void h5load(
    const string &file, 
    const string &dataset
  )
  {
    notice_macro("about to open file: " << file)
    h5f.openFile(file, H5F_ACC_RDONLY);

    notice_macro("about to read dataset: " << dataset)
    h5d = h5f.openDataSet(dataset);
    h5s = h5d.getSpace();
  }

  public:

  //ctor
  PlotterCommon(const string &file):
    file(file)
  {
    // init dt and outfreq
    {
      notice_macro("about to open file: " << file << "/const.h5")
      h5f.openFile(file + "/const.h5", H5F_ACC_RDONLY);

      notice_macro("about to read group: advection")
      h5g = h5f.openGroup("advection");

      float dt;
      {
        auto attr = h5g.openAttribute("dt");
        attr.read(attr.getDataType(), &dt);
      }
      map["dt"] = dt;


      // read number of timesteps
      hsize_t n;
      h5load(file + "/const.h5", "T");
      h5s.getSimpleExtentDims(&n, NULL);
      this->map["t"] = n;
      // read timesteps
      timesteps.resize(n);
      h5d.read(timesteps.data(), H5::PredType::NATIVE_FLOAT);
      // read output frequency
      float outfreq;
      {
        auto root_group = h5f.openGroup("/");
        auto attr = root_group.openAttribute("user_params outfreq");
        attr.read(attr.getDataType(), &outfreq);
      }
      map["outfreq"] = outfreq;
    }
  }
};

