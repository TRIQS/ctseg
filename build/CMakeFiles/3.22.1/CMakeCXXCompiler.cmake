set(CMAKE_CXX_COMPILER "/cm/shared/sw/nix/store/i3n3fdm2hkrn8njaklii8ni1g9lp7imz-llvm-13.0.0/bin/clang++")
set(CMAKE_CXX_COMPILER_ARG1 "")
set(CMAKE_CXX_COMPILER_ID "Clang")
set(CMAKE_CXX_COMPILER_VERSION "13.0.0")
set(CMAKE_CXX_COMPILER_VERSION_INTERNAL "")
set(CMAKE_CXX_COMPILER_WRAPPER "")
set(CMAKE_CXX_STANDARD_COMPUTED_DEFAULT "14")
set(CMAKE_CXX_EXTENSIONS_COMPUTED_DEFAULT "ON")
set(CMAKE_CXX_COMPILE_FEATURES "cxx_std_98;cxx_template_template_parameters;cxx_std_11;cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates;cxx_std_14;cxx_aggregate_default_initializers;cxx_attribute_deprecated;cxx_binary_literals;cxx_contextual_conversions;cxx_decltype_auto;cxx_digit_separators;cxx_generic_lambdas;cxx_lambda_init_captures;cxx_relaxed_constexpr;cxx_return_type_deduction;cxx_variable_templates;cxx_std_17;cxx_std_20;cxx_std_23")
set(CMAKE_CXX98_COMPILE_FEATURES "cxx_std_98;cxx_template_template_parameters")
set(CMAKE_CXX11_COMPILE_FEATURES "cxx_std_11;cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates")
set(CMAKE_CXX14_COMPILE_FEATURES "cxx_std_14;cxx_aggregate_default_initializers;cxx_attribute_deprecated;cxx_binary_literals;cxx_contextual_conversions;cxx_decltype_auto;cxx_digit_separators;cxx_generic_lambdas;cxx_lambda_init_captures;cxx_relaxed_constexpr;cxx_return_type_deduction;cxx_variable_templates")
set(CMAKE_CXX17_COMPILE_FEATURES "cxx_std_17")
set(CMAKE_CXX20_COMPILE_FEATURES "cxx_std_20")
set(CMAKE_CXX23_COMPILE_FEATURES "cxx_std_23")

set(CMAKE_CXX_PLATFORM_ID "Linux")
set(CMAKE_CXX_SIMULATE_ID "")
set(CMAKE_CXX_COMPILER_FRONTEND_VARIANT "GNU")
set(CMAKE_CXX_SIMULATE_VERSION "")




set(CMAKE_AR "/cm/shared/sw/nix/store/i3n3fdm2hkrn8njaklii8ni1g9lp7imz-llvm-13.0.0/bin/llvm-ar")
set(CMAKE_CXX_COMPILER_AR "/cm/shared/sw/nix/store/i3n3fdm2hkrn8njaklii8ni1g9lp7imz-llvm-13.0.0/bin/llvm-ar")
set(CMAKE_RANLIB "/cm/shared/sw/nix/store/i3n3fdm2hkrn8njaklii8ni1g9lp7imz-llvm-13.0.0/bin/llvm-ranlib")
set(CMAKE_CXX_COMPILER_RANLIB "/cm/shared/sw/nix/store/i3n3fdm2hkrn8njaklii8ni1g9lp7imz-llvm-13.0.0/bin/llvm-ranlib")
set(CMAKE_LINKER "/cm/shared/sw/nix/store/i3n3fdm2hkrn8njaklii8ni1g9lp7imz-llvm-13.0.0/bin/ld.lld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCXX )
set(CMAKE_CXX_COMPILER_LOADED 1)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_CXX_ABI_COMPILED TRUE)

set(CMAKE_CXX_COMPILER_ENV_VAR "CXX")

set(CMAKE_CXX_COMPILER_ID_RUN 1)
set(CMAKE_CXX_SOURCE_FILE_EXTENSIONS C;M;c++;cc;cpp;cxx;m;mm;mpp;CPP;ixx;cppm)
set(CMAKE_CXX_IGNORE_EXTENSIONS inl;h;hpp;HPP;H;o;O;obj;OBJ;def;DEF;rc;RC)

foreach (lang C OBJC OBJCXX)
  if (CMAKE_${lang}_COMPILER_ID_RUN)
    foreach(extension IN LISTS CMAKE_${lang}_SOURCE_FILE_EXTENSIONS)
      list(REMOVE_ITEM CMAKE_CXX_SOURCE_FILE_EXTENSIONS ${extension})
    endforeach()
  endif()
endforeach()

set(CMAKE_CXX_LINKER_PREFERENCE 30)
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 1)

# Save compiler ABI information.
set(CMAKE_CXX_SIZEOF_DATA_PTR "8")
set(CMAKE_CXX_COMPILER_ABI "ELF")
set(CMAKE_CXX_BYTE_ORDER "LITTLE_ENDIAN")
set(CMAKE_CXX_LIBRARY_ARCHITECTURE "")

if(CMAKE_CXX_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CXX_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CXX_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CXX_COMPILER_ABI}")
endif()

if(CMAKE_CXX_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_CXX_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES "/cm/shared/sw/nix/store/hryq9drhdswvw0khbqaqk9p5k20vpzs3-intel-oneapi-mkl-2022.0.1/mkl/2022.0.1/include;/cm/shared/sw/nix/store/clmmc53sav2i4hvg53q5866nd9nygc4b-python-3.9.9-view/include/python3.9;/cm/shared/sw/nix/store/i3n3fdm2hkrn8njaklii8ni1g9lp7imz-llvm-13.0.0/include;/cm/shared/sw/nix/store/clmmc53sav2i4hvg53q5866nd9nygc4b-python-3.9.9-view/include;/cm/shared/sw/nix/store/v92lzmn32qlq3l7hrah0rzkazpcr7m8k-boost-1.78.0/include;/cm/shared/sw/nix/store/2nq01a1xf3r3lgjbgp8wp9xpx2y8hcj5-hdf5-1.10.8/include;/cm/shared/sw/nix/store/wf2b4c91f1s7n4cysvk0bshfnywvasmv-nfft-3.5.2/include;/cm/shared/sw/nix/store/avq18lr6jgg6fpikan0cdzaf61xzigkk-fftw-3.3.10/include;/cm/shared/sw/nix/store/j4skrrvqyvk9bc2p7zqjd0qsfg4sg8nd-gmp-6.2.1/include;/cm/shared/sw/nix/store/dvc90n2l1p8lz5an0d4lnk2ahca50hfm-openmpi-4.0.7/include;/cm/shared/sw/nix/store/ddhcca3mhxdlqwm95fklqgv3ximkkwgz-flexiblas-3.0.4/include;/cm/shared/sw/nix/store/aqkirva9k1vggmfcv7hy9rsy0kavh56m-gcc-10.3.0/include;/cm/shared/apps/slurm/current/include;/cm/shared/sw/nix/store/i3n3fdm2hkrn8njaklii8ni1g9lp7imz-llvm-13.0.0/include/c++/v1;/cm/shared/sw/nix/store/i3n3fdm2hkrn8njaklii8ni1g9lp7imz-llvm-13.0.0/lib/clang/13.0.0/include;/usr/local/include;/usr/include")
set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "c++;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/cm/shared/sw/nix/store/hmphldjnv7440cwpv29ph5d5bp3lbxz5-gcc-11.2.0/lib/gcc/x86_64-pc-linux-gnu/11.2.0;/cm/shared/sw/nix/store/hmphldjnv7440cwpv29ph5d5bp3lbxz5-gcc-11.2.0/lib64;/lib64;/usr/lib64;/cm/shared/sw/nix/store/i3n3fdm2hkrn8njaklii8ni1g9lp7imz-llvm-13.0.0/lib;/lib;/usr/lib;/cm/shared/sw/nix/store/hryq9drhdswvw0khbqaqk9p5k20vpzs3-intel-oneapi-mkl-2022.0.1/mkl/2022.0.1/lib/intel64;/cm/shared/sw/nix/store/hryq9drhdswvw0khbqaqk9p5k20vpzs3-intel-oneapi-mkl-2022.0.1/lib;/cm/shared/sw/nix/store/bcx7ls8qbd098qhkg9axnpdpsavb287l-py-mpi4py-3.1.2-view/lib;/cm/shared/sw/nix/store/clmmc53sav2i4hvg53q5866nd9nygc4b-python-3.9.9-view/lib;/cm/shared/sw/nix/store/v92lzmn32qlq3l7hrah0rzkazpcr7m8k-boost-1.78.0/lib64;/cm/shared/sw/nix/store/v92lzmn32qlq3l7hrah0rzkazpcr7m8k-boost-1.78.0/lib;/cm/shared/sw/nix/store/2nq01a1xf3r3lgjbgp8wp9xpx2y8hcj5-hdf5-1.10.8/lib64;/cm/shared/sw/nix/store/2nq01a1xf3r3lgjbgp8wp9xpx2y8hcj5-hdf5-1.10.8/lib;/cm/shared/sw/nix/store/wf2b4c91f1s7n4cysvk0bshfnywvasmv-nfft-3.5.2/lib;/cm/shared/sw/nix/store/avq18lr6jgg6fpikan0cdzaf61xzigkk-fftw-3.3.10/lib64;/cm/shared/sw/nix/store/avq18lr6jgg6fpikan0cdzaf61xzigkk-fftw-3.3.10/lib;/cm/shared/sw/nix/store/j4skrrvqyvk9bc2p7zqjd0qsfg4sg8nd-gmp-6.2.1/lib;/cm/shared/sw/nix/store/dvc90n2l1p8lz5an0d4lnk2ahca50hfm-openmpi-4.0.7/lib;/cm/shared/sw/nix/store/ddhcca3mhxdlqwm95fklqgv3ximkkwgz-flexiblas-3.0.4/lib64;/cm/shared/sw/nix/store/aqkirva9k1vggmfcv7hy9rsy0kavh56m-gcc-10.3.0/lib64;/cm/shared/sw/nix/store/aqkirva9k1vggmfcv7hy9rsy0kavh56m-gcc-10.3.0/lib;/cm/shared/apps/slurm/current/lib64")
set(CMAKE_CXX_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
