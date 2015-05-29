 
	function(GetLastPart path name)
		string(FIND ${path} "/" last_slash REVERSE)
		string(LENGTH ${path} dir_length)
		MATH(EXPR sub_str_length "${dir_length} - ${last_slash} - 1" )
		MATH(EXPR last_slash "${last_slash} + 1")
		string(SUBSTRING ${path} ${last_slash} ${sub_str_length} ${name})
		set(${name} ${${name}} PARENT_SCOPE)
	endfunction()

	file(GLOB dirs *)
	set(local_inc_dirs "")
	foreach(dir ${dirs})
	  if (IS_DIRECTORY ${dir})
	    include_directories(${dir})
	    set(local_inc_dirs ${local_inc_dirs} ${dir})
	  endif()
	endforeach(dir)
	set(INC_DIRS ${INC_DIRS} ${local_inc_dirs} PARENT_SCOPE)
	
set(DIRS
  pbd_rod
  pose
  subspace
  volumetric_rod
  tet_mesh_io
  tri_mesh_io
  vega
  window
  global
  qt_opengl_util
  bvh
  lib
)

	set(local_libs "")
	foreach(dir ${DIRS})
	  set(local_libs ${local_libs} ${dir})
	  add_subdirectory(${dir})
	endforeach(dir)
	set(LIBS ${LIBS} ${local_libs} PARENT_SCOPE)
	
