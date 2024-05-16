add_library(addimages SHARED "addimages.cpp")
add_library(joinrgb   SHARED "joinrgb.cpp"  )
add_library(cylclip   SHARED "cylclip.cpp"  )
add_library(kdbench   SHARED "kdbench.cpp"  )
add_library(tonemap   SHARED "tonemap.cpp"  )

add_library(matgen   SHARED "mat_gen.cpp"  )

target_link_libraries(addimages mitsuba-python)
target_link_libraries(joinrgb mitsuba-python)
target_link_libraries(cylclip mitsuba-python)
target_link_libraries(kdbench getopt mitsuba-python)
target_link_libraries(tonemap getopt mitsuba-python)

target_link_libraries(matgen mitsuba-python)

set_target_properties(addimages PROPERTIES FOLDER "plugins")
set_target_properties(joinrgb PROPERTIES FOLDER "plugins")
set_target_properties(cylclip PROPERTIES FOLDER "plugins")
set_target_properties(kdbench PROPERTIES FOLDER "plugins")
set_target_properties(tonemap PROPERTIES FOLDER "plugins")

set_target_properties(matgen PROPERTIES FOLDER "plugins")