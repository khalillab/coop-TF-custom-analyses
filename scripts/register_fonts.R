
main = function(fonts_path,
                output_path){
    library(extrafont)

    ttf_import(fonts_path)
    loadfonts()

    file.create(output_path)
}

main(fonts_path = snakemake@input[["fonts_path"]],
     output_path = snakemake@output[["output_path"]])

