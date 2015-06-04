# UCSCBasics.py
#
# Description:
#  Hooks into the UCSC database
#

class URLfactory:
  # init
  # Pre: db database ie hg19
  def __init__(self,db):
    self.db = db
    self.session = False

  def set_session(self,session_string):
    self.session = session_string

  # for the base-0 start base-1 stop type
  def url_from_bed_coordinates(self,chr,start0,stop1):
    #https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr12%3A6643585-6647537
    url =  'https://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.db+'&'
    url += 'position='+chr+'%3A'+str(start0+1)+'-'+str(stop1)
    if self.session:
      url += '&hgsid='+self.session
    return url

  def excel_from_bed_coordinates(self,chr,start0,stop1):
    field =  '=HYPERLINK("'+self.url_from_bed_coordinates(chr,start0,stop1)+'",'
    field += '"'+chr+':'+str(start0+1)+'-'+str(stop1)+'")'
    return field
