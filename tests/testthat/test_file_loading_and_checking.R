context("")

write.table(DNABarcodeCompatibility::IlluminaIndexesRaw ,
            textfile <- tempfile(),
            row.names = FALSE, col.names = FALSE, quote=FALSE)

test_that("file_loading_and_checking loads file", {
    expect_is(file_loading_and_checking(textfile), "data.frame")
})
