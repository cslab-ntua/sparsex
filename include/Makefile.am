#
# Makefile.am -- SparseX library
#
# Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA
# Copyright (C) 2014, Vasileios Karakasis
# All rights reserved.
#
# This file is distributed under the BSD License. See LICENSE.txt for details.
#

RM = /bin/rm -f
install_data = @INSTALL_DATA@

nobase_include_HEADERS = \
	sparsex/common.h \
	sparsex/error.h \
	sparsex/matvec.h \
	sparsex/sparsex.h \
	sparsex/timing.h \
	sparsex/types.h \
	sparsex/internals/cdecl.h \
	sparsex/internals/futex.h \
	sparsex/internals/numa_util.h \
	sparsex/internals/Affinity.hpp \
	sparsex/internals/Allocators.hpp \
	sparsex/internals/Barrier.hpp \
	sparsex/internals/Csr.hpp \
	sparsex/internals/CodeExecutor.hpp \
	sparsex/internals/Csx.hpp \
	sparsex/internals/CsxBuild.hpp \
	sparsex/internals/CsxGetSet.hpp \
	sparsex/internals/CsxJit.hpp \
	sparsex/internals/CsxKernels.hpp \
	sparsex/internals/CsxManager.hpp \
	sparsex/internals/CsxSaveRestore.hpp \
	sparsex/internals/CsxSpmv.hpp \
	sparsex/internals/CsxUtil.hpp \
	sparsex/internals/CtlBuilder.hpp \
	sparsex/internals/CtlUtil.hpp \
	sparsex/internals/Delta.hpp \
	sparsex/internals/DynamicArray.hpp \
	sparsex/internals/Element.hpp \
	sparsex/internals/EncodingManager.hpp \
	sparsex/internals/Encodings.hpp \
	sparsex/internals/Facade.hpp \
	sparsex/internals/JitCompiler.hpp \
	sparsex/internals/JitConfig.hpp \
	sparsex/internals/JitUtil.hpp \
	sparsex/internals/Map.hpp \
	sparsex/internals/Mmf.hpp \
	sparsex/internals/Node.hpp \
	sparsex/internals/Rcm.hpp \
	sparsex/internals/Runtime.hpp \
	sparsex/internals/SparseInternal.hpp \
	sparsex/internals/SparseMatrix.hpp \
	sparsex/internals/SparsePartition.hpp \
	sparsex/internals/SpmMt.hpp \
	sparsex/internals/SpmvMethod.hpp \
	sparsex/internals/Statistics.hpp \
	sparsex/internals/TemplateText.hpp \
	sparsex/internals/ThreadPool.hpp \
	sparsex/internals/Timer.hpp \
	sparsex/internals/TimerCollection.hpp \
	sparsex/internals/Types.hpp \
	sparsex/internals/Utility.hpp \
	sparsex/internals/Vector.hpp \
	sparsex/internals/Xform.hpp \
	sparsex/internals/logger/Logger.hpp \
	sparsex/internals/logger/LoggerUtil.hpp \
	sparsex/internals/logger/OutputPolicy.hpp

# Fill in the correct paths before installing
install-data-local:
	$(SED) -e "/.*#define SPX_JIT_INCLUDE/ s/\".*\"/\"`echo $(includedir) | $(SED) -e 's/\//\\\\\//g'`\"/g; /.*#define SPX_MULT_TEMPLATE_DIR/ s/\".*\"/\"`echo $(datarootdir) | $(SED) -e 's/\//\\\\\//g'`\/sparsex\/templates\"/g" sparsex/internals/Config.hpp > sparsex/internals/Config.hpp.install; \
test -e $(includedir)/sparsex/internals || \
$(MKDIR_P) $(includedir)/sparsex/internals; \
$(install_data) sparsex/internals/Config.hpp.install $(includedir)/sparsex/internals/Config.hpp; \
$(install_data) sparsex/config.h $(includedir)/sparsex/config.h

uninstall-local:
	$(RM) $(includedir)/sparsex/config.h
	$(RM) $(includedir)/sparsex/internals/Config.hpp
