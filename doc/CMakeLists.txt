add_custom_command(OUTPUT userguide.html
                   COMMAND pandoc -r markdown -t html -o userguide.html -c ${CMAKE_CURRENT_SOURCE_DIR}/parosol.css -S -N -s --toc ${CMAKE_CURRENT_SOURCE_DIR}/userguide.mkd DEPENDS userguide.mkd parosol.css
                   )

add_custom_target(userguide ALL DEPENDS userguide.html)

install(FILES userguide.html DESTINATION share/doc/parosol)
install(FILES parosol.css DESTINATION share/doc/parosol)
