==================
1. Getting started
==================

------------
Installation
------------

^^^^^^^^^^^^
Dependencies
^^^^^^^^^^^^

^^^^^^^^^^^^^^^^
Download & Build
^^^^^^^^^^^^^^^^

-------------------------
Using *obscura* as a tool
-------------------------

^^^^^^^^^^^^^^^^^^^^^^
The configuration file
^^^^^^^^^^^^^^^^^^^^^^


----------------------------
Using *obscura* as a library
----------------------------

.. ## Getting started

.. <details><summary>1. Dependencies</summary>

.. <details><summary>[boost](https://www.boost.org/)</summary>
.. To install `boost` on a Mac using [homebrew](https://brew.sh/), simply run
.. ```
.. >brew install boost
.. ```

.. On Linux machines, run
.. ```
.. >sudo apt-get update && sudo apt-get install -yq libboost-all-dev
.. ```
.. </p>
.. </details>

.. <details><summary>[libconfig](https://hyperrealm.github.io/libconfig/)</summary>
.. To install `libconfig` on a Mac using [homebrew](https://brew.sh/), simply run
.. ```
.. >brew install libconfig
.. ```

.. On Linux machines, you can build `libconfig` via
.. ```
.. >wget https://hyperrealm.github.io/libconfig/dist/libconfig-1.7.2.tar.gz
.. >tar -xvzf libconfig-1.7.2.tar.gz
.. >pushd libconfig-1.7.2
.. >./configure
.. >make
.. >sudo make install
.. >popd
.. ```

.. </p>
.. </details>

.. <details><summary>[libphysica](https://github.com/temken/libphysica)</summary>
.. `libphysica` does not need to be installed. It will be downloaded and compiled during the CMake build.
.. </p>
.. </details>

.. </p>
.. </details>

.. <details><summary>2. Download & Installation</summary>
.. The `obscura` source code can be downloaded by cloning this git repository:

.. ```
.. >git clone https://github.com/temken/`obscura`.git 
.. >cd obscura
.. ```

.. The code is compiled and the executable is created using CMake.

.. ```
.. >cmake -E make_directory build
.. >cd build
.. >cmake -DCMAKE_BUILD_TYPE=Release -DCODE_COVERAGE=OFF ..
.. >cmake --build . --config Release
.. >cmake --install .
.. ```

.. If everything worked well, there should be the executable *obscura* in the */bin/* folder.

.. </p>
.. </details>

.. <details><summary>3. Usage as a tool</summary>

.. </p>
.. </details>

.. <details><summary>4. Usage as an external library</summary>

.. </p>
.. </details>