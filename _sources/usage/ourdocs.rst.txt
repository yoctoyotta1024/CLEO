Our Docs
========

.. note::
   Please consider that this project is under active development
   and our documentation is a work in progress.

We use Sphinx to build our documentation. Through the Breathe
extension we integrate C++ documentation created using Doxygen. You can
build our documentation locally by building the Doxygen .xml files
followed by the Sphinx .html files, which you can then view
in your preferred browser. E.g.

.. code-block:: console

  $ cd ~/CLEO/docs && mkdir build && mkdir build/doxygen
  $ doxygen doxygen/doxygen.dox && make html
  $ open build/html/index.html


If you would like to contribute to the documentation, have questions or
would like clarification or addition to the documentation about
a particular topic, please :ref:`contact us! <contact>` or `open a new
issue <https://github.com/yoctoyotta1024/CLEO/issues/new>`_ on
our GitHub repository.
