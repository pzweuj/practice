# coding=utf-8
# pzw
# 20190806

from PIL import Image, ImageDraw


# def add_watermark_to_image(image, watermark):
# 	rgba_image = image.convert("RGBA")
# 	rgba_watermark = Image.open(watermark).convert("RGBA")
 
# 	image_x, image_y = rgba_image.size
# 	watermark_x, watermark_y = rgba_watermark.size
 
# 	# # 缩放图片
# 	# scale = 10
# 	# watermark_scale = max(image_x / (scale * watermark_x), image_y / (scale * watermark_y))
# 	# new_size = (int(watermark_x * watermark_scale), int(watermark_y * watermark_scale))
# 	# rgba_watermark = rgba_watermark.resize(new_size, resample=Image.ANTIALIAS)
# 	# 透明度
# 	rgba_watermark_mask = rgba_watermark.convert("L").point(lambda x: 90)
# 	rgba_watermark.putalpha(rgba_watermark_mask)
 
# 	watermark_x, watermark_y = rgba_watermark.size
# 	# 水印位置
# 	rgba_image.paste(rgba_watermark, (2500, 2000), rgba_watermark_mask)



img = Image.open("image41.jpeg")
a = img.size

cut1 = img.crop((0, 100, a[0]/2-300, a[1]))
# cut1.save("test.jpeg")

cut2 = img.crop((a[0]/2+220, 100, a[0], a[1]))
# cut2_draw = ImageDraw.Draw(cut2)

mark = Image.open("fix1.jpg")

layer = Image.new("RGBA", cut2.size, (0, 0, 0, 0))
layer.paste(mark, (1600, 3800))
out = Image.composite(layer, cut2, layer)
# out.show()
out.save("test2.jpeg")

layer1 = Image.new("RGBA", cut1.size, (0, 0, 0, 0))
layer1.paste(mark, (1500, 3600))
out1 = Image.composite(layer1, cut1, layer1)
# out1.show()
out1.save("test1.jpeg")