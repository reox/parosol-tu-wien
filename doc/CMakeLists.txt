set (OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})



add_custom_command(OUTPUT userguide.html
                   COMMAND pandoc
                   ARGS -r markdown -t html -o userguide.html -c parosol.css -S -N -s --toc userguide.mkd
                   MAIN_DEPENDENCY userguide.mkd
                   DEPENDS parosol.css
                   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                   )

               #install (FILES userguide.html DESTINATION ${TARGET_DOCUMENTATION_TEXT_FOLDER})
