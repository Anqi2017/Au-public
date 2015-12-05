import zlib, json, base64

def encode_32(obj):
  compressed_string = zlib.compress(json.dumps(obj),9)
  enc_string = base64.b32encode(compressed_string)
  return 'SER32_'+enc_string.rstrip('=')

def decode_32(safename):
  frag = safename.lstrip('SER32_')
  padding = (8-(len(frag) % 8)) % 8
  c = base64.b32decode(frag+'='*padding)
  return json.loads(zlib.decompress(c))

def encode_64(obj):
  compressed_string = zlib.compress(json.dumps(obj),9)
  enc_string = base64.b64encode(compressed_string)
  return 'SER64_'+enc_string

def decode_64(safename):
  frag = safename.lstrip('SER64_')
  c = base64.b64decode(frag)
  return json.loads(zlib.decompress(c))
