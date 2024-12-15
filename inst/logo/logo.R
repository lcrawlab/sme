library(hexSticker)
library(magick)
imgurl <- image_read("favicon.png")
sysfonts::font_add_google("Zilla Slab", "pf", regular.wt = 500)

r <- hexSticker::sticker(
  subplot = imgurl,
  s_x = 1,
  s_y = 1.25,
  s_width = 1.25,
  s_height = 1.25,
  package = "mmeR",
  p_x = 1,
  p_y = 0.6,
  p_size = 25,
  h_size = 1.2,
  p_family = "pf",
  p_color = "black",
  h_fill = "white", #"#FFF9F2",
  h_color = "#7562A4",
  dpi = 320,
  filename = "logo.png"
)

plot(r)
