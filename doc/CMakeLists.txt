# FIXME
add_custom_command(OUTPUT userguide.html
                   COMMAND pandoc -r markdown -t html -o userguide.html -c parosol.css -S -N -s --toc userguide.mkd
                   )

