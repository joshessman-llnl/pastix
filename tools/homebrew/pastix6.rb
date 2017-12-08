class Pastix6 < Formula
  desc "Parallel solver for sparse linear systems based on direct methods"
  homepage "https://gitlab.inria.fr/solverstack/pastix"
  url "https://www.labri.fr/perso/ramet/pastix-master-6.0.0.tar.bz2"
  sha256 "fe365c01e00eb4918ff8dd2d9720df500d461802179a625afce405a0e600c42f"

  head "https://gitlab.inria.fr/solverstack/pastix.git"

  bottle :disable, "needs to be rebuilt"

  depends_on "openblas"
  depends_on "hwloc"              # Could be optinal but strongly recommanded
  depends_on "scotch"             # Could be optinal but strongly recommanded
  depends_on "metis" => :optional # Use METIS ordering.
  depends_on "gcc"   => :build    # GNU Fortran is now provided as part of GCC
  depends_on "cmake" => :build

  conflicts_with "pastix", :because => "conflicts with default links"

  def install
    args = ["-DCMAKE_INSTALL_PREFIX=#{prefix}",
            "-DBUILD_SHARED_LIBS=ON",
            "-DBUILD_DOCUMENTATION=OFF",
            "-DBUILD_64bits=OFF",
            "-DPASTIX_INT64=OFF",
            "-DPASTIX_ORDERING_SCOTCH=ON",
            "-DPASTIX_WITH_FORTRAN=ON",
            "-DPASTIX_WITH_MPI=OFF",
            "-DPASTIX_WITH_CUDA=OFF",
            "-DPASTIX_WITH_STARPU=OFF",
            "-DPASTIX_WITH_PARSEC=OFF"]
    args += ["-DPASTIX_ORDERING_METIS=ON"] if build.with? "metis"
    mkdir "build" do
      system "cmake", "..", *args
      system "make"
      system "make", "install"
    end
    pkgshare.install "example" # Contains all test programs.
  end

  def caveats; <<-EOS.undent
    Set the PYTHONPATH environment variable:
      export PYTHONPATH=#{prefix}/lib/python:$PYTHONPATH
    Try python example with:
      python #{prefix}/examples/simple.py
    Or simple example with:
      #{prefix}/examples/simple -9 10:10:10
    EOS
  end

  test do
    system "#{prefix}/examples/simple", "-9", "10:10:10"
    ohai "All test output is in #{HOMEBREW_LOGS}/pastix. Please check."
  end
end
