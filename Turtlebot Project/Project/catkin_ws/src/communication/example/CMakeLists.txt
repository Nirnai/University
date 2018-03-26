set(MOC_H
  testgui.h
)

set(SRC_H
  testgui.h
)


set(SRC_CPP
  testgui.cpp
  main.cpp
)

QT4_WRAP_CPP( MOC_CPP ${MOC_H} )

add_executable( testgui
	${SRC_CPP}
  ${SRC_H} 
	${MOC_CPP}
)
  
target_link_libraries( testgui         
	${QT_LIBRARIES}
  referee 
) 
