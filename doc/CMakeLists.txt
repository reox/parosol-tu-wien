# FIXME
add_custom_command(OUTPUT userguide.html
                   COMMAND pandoc -r markdown -t html -o userguide.html -c parosol.css -S -N -s --toc userguide.mkd)

add_custom_target(userguide ALL DEPENDS userguide.html)

install(FILES userguide.html DESTINATION share/doc/parosol)

