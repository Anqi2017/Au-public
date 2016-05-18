#!/usr/bin/python
import sys, argparse, base64, re, os


def main():
  parser = argparse.ArgumentParser(description="Put css style sheets and PNG images into html file",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Input file")
  args = parser.parse_args()
  
  filecontents = open(args.input).read()
  imgs = []
  for m in re.finditer('(<\s*.*img.*>)',filecontents):
    imgs.append([m.start(),m.end()])
  prev = 0
  newcontents = ''
  for i in range(len(imgs)):
    newcontents += filecontents[prev:imgs[i][0]]
    newcontents += do_image_tag(filecontents[imgs[i][0]:imgs[i][1]],args)
    prev = imgs[i][1]
  newcontents += filecontents[prev:]
  styles = []
  for m in re.finditer('(<\s*.*type.*text/css.*>)',newcontents):
    styles.append([m.start(),m.end()])
  prev = 0    
  for i in range(len(styles)):
    print newcontents[prev:styles[i][0]]
    print do_style_sheet(newcontents[styles[i][0]:styles[i][1]],args)
    prev = styles[i][1]
  print newcontents[prev:]


def do_style_sheet(style_sheet,args):
  m=re.match('^(.*)(href\s*=\s*["\'][^"\']*["\'])(.*)$',style_sheet)
  if not m: 
    return style_sheet #cant replace for some reason
  start = m.group(1)
  finish = m.group(3)
  src_full = m.group(2)
  m = re.match('href\s*=\s*["\']([^"\']+)["\']',src_full)
  if not m:
    return style_sheet #cant replace for some reason
  srcpathpart = m.group(1)
  srcpath = os.path.dirname(args.input)+'/'+srcpathpart
  if not re.search('\.css',srcpath): return style_sheet
  if not os.path.isfile(srcpath): return style_sheet
  encoded = base64.b64encode(open(srcpath,'rb').read())
  disabler = "a {\n  pointer-events: none;\n}\n"
  return '<style>'+disabler+"\n"+open(srcpath,'rb').read()+'</style>'

def do_image_tag(img_tag,args):
  m=re.match('^(.*)(src\s*=\s*["\'][^"\']*["\'])(.*)$',img_tag)
  if not m: 
    return img_tag #cant replace for some reason
  start = m.group(1)
  finish = m.group(3)
  src_full = m.group(2)
  m = re.match('src\s*=\s*["\']([^"\']+)["\']',src_full)
  if not m:
    return img_tag #cant replace for some reason
  srcpathpart = m.group(1)
  srcpath = os.path.dirname(args.input)+'/'+srcpathpart
  if not re.search('\.png',srcpath): return img_tag
  if not os.path.isfile(srcpath): return img_tag
  encoded = base64.b64encode(open(srcpath,'rb').read())
  return start+' src="data:image/png;base64,'+encoded+'" '+finish

if __name__=="__main__":
  main()
