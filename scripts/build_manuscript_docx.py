"""Build the editable manuscript from Markdown and figure files."""
from pathlib import Path
import re,argparse
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
doc=Document();sec=doc.sections[0];sec.page_width=Inches(8.5);sec.page_height=Inches(11);sec.top_margin=sec.bottom_margin=Inches(.8);sec.left_margin=sec.right_margin=Inches(.85)
for name in ['Normal','Title','Heading 1','Heading 2','Heading 3','Caption']:
 st=doc.styles[name];st.font.name='Times New Roman';st.font.color.rgb=RGBColor(0,0,0);st.font.size=Pt(11)
 st.paragraph_format.space_after=Pt(6);st.paragraph_format.line_spacing=1.1
for name,size in [('Title',18),('Heading 1',14),('Heading 2',12),('Heading 3',11)]:
 st=doc.styles[name];st.font.size=Pt(size);st.font.bold=True;st.paragraph_format.keep_with_next=True;st.paragraph_format.space_before=Pt(10)
doc.styles['Title'].paragraph_format.space_before=Pt(0)
doc.styles['Caption'].font.size=Pt(10)
doc.styles['Caption'].font.bold=False
for el in list(doc.styles.element.iter(qn('w:pBdr'))):el.getparent().remove(el)
for st in doc.styles:
 for el in st.element.iter(qn('w:rFonts')):
  for attr in list(el.attrib):
   if attr.endswith('Theme'):del el.attrib[attr]
doc.core_properties.title='PocketToDynamics technical manuscript';doc.core_properties.author='Juan Philippe Teixeira; Daniel Ferreira de Lima Neto; Miklos Maximiliano Bajay'
# Page numbering aids scientific review without a decorative header.
pf=sec.footer.paragraphs[0];pf.alignment=WD_ALIGN_PARAGRAPH.CENTER
fld=OxmlElement('w:fldSimple');fld.set(qn('w:instr'),'PAGE');pf._p.append(fld)
def inline(p,text):
 for token in re.split(r'(\*\*.*?\*\*|\*[^*]+\*|`[^`]+`|\[[^\]]+\]\((?:[^()]|\([^()]*\))+\))',text):
  if token.startswith('**') and token.endswith('**'):r=p.add_run(token[2:-2]);r.bold=True
  elif token.startswith('`') and token.endswith('`'):r=p.add_run(token[1:-1]);r.font.name='Consolas';r.font.size=Pt(9)
  elif re.fullmatch(r'\[[^\]]+\]\((?:[^()]|\([^()]*\))+\)',token):
   m=re.match(r'\[([^\]]+)\]\(((?:[^()]|\([^()]*\))+)\)',token)
   rel=p.part.relate_to(m[2],'http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink',is_external=True);link=OxmlElement('w:hyperlink');link.set(qn('r:id'),rel);rr=OxmlElement('w:r');tt=OxmlElement('w:t');tt.text=m[1];rr.append(tt);link.append(rr);p._p.append(link)
  elif token.startswith('*') and token.endswith('*'):r=p.add_run(token[1:-1]);r.italic=True
  else:p.add_run(token)
def table(lines):
 rows=[[c.strip() for c in l.strip().strip('|').split('|')] for l in lines];rows=[r for r in rows if not all(re.fullmatch(r':?-+:?',x) for x in r)]
 t=doc.add_table(rows=0,cols=len(rows[0]));t.alignment=WD_TABLE_ALIGNMENT.CENTER;t.autofit=False
 widths=([1.0,2.65,3.15] if len(rows[0])==3 else [6.8/len(rows[0])]*len(rows[0]))
 for col,w in zip(t.columns,widths):col.width=Inches(w)
 for i,row in enumerate(rows):
  cells=t.add_row().cells
  trpr=t.rows[-1]._tr.get_or_add_trPr();nosplit=OxmlElement('w:cantSplit');trpr.append(nosplit)
  if i==0:el=OxmlElement('w:tblHeader');trpr.append(el)
  for c,txt,w in zip(cells,row,widths):
   c.width=Inches(w);c.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER;pr=c._tc.get_or_add_tcPr();b=OxmlElement('w:tcBorders')
   for side in ['top','bottom']:
    e=OxmlElement('w:'+side);e.set(qn('w:val'),'single');e.set(qn('w:sz'),'6' if i==0 else '3');e.set(qn('w:color'),'7F7F7F' if i==0 else 'D9D9D9');b.append(e)
   for side in ['left','right','insideV']:
    e=OxmlElement('w:'+side);e.set(qn('w:val'),'nil');b.append(e)
   pr.append(b)
   if i==0:sh=OxmlElement('w:shd');sh.set(qn('w:fill'),'E6E6E6');pr.append(sh)
   pp=c.paragraphs[0];pp.paragraph_format.space_after=Pt(4);pp.paragraph_format.space_before=Pt(4);inline(pp,txt)
   for r in pp.runs:r.font.size=Pt(9.5);r.bold=i==0
 doc.add_paragraph().paragraph_format.space_after=Pt(0)
lines=a.source.read_text().splitlines();i=0;pending_image=None;in_references=False
while i<len(lines):
 l=lines[i];i+=1
 if not l.strip():continue
 if l.startswith('!['):
  m=re.match(r'!\[([^\]]*)\]\(((?:[^()]|\([^()]*\))+)\)',l);pending_image=(m[1],a.source.parent/m[2]);continue
 if l.startswith('*Figure '):
  caption=l.strip('*')
  if pending_image is not None:
   alt,path=pending_image;assert path.exists(),path
   pp=doc.add_paragraph();pp.paragraph_format.keep_with_next=True;pp.paragraph_format.space_before=Pt(4);pp.paragraph_format.space_after=Pt(2);pp.add_run().add_picture(str(path.resolve()),width=Inches(6.55))
   for node in pp._p.iter(qn('wp:docPr')):node.set('descr',alt)
   cp=doc.add_paragraph(caption,style='Caption');cp.paragraph_format.keep_together=True;cp.paragraph_format.space_after=Pt(9)
  pending_image=None;continue
 if l.startswith('```'):
  lang=l[3:];code=[]
  while i<len(lines) and not lines[i].startswith('```'):code.append(lines[i]);i+=1
  i+=1
  if lang=='mermaid':continue
  pp=doc.add_paragraph();pp.paragraph_format.keep_together=True
  for j,ln in enumerate(code):r=pp.add_run(ln+('\n' if j<len(code)-1 else ''));r.font.name='Consolas';r.font.size=Pt(9)
  continue
 if l=='\\[':
  while i<len(lines) and lines[i]!='\\]':i+=1
  i+=1;pp=doc.add_paragraph();pp.alignment=WD_ALIGN_PARAGRAPH.CENTER;math=OxmlElement('m:oMath');rr=OxmlElement('m:r');tt=OxmlElement('m:t');tt.text='Fᵢ = −RT ln(Pᵢ/Pₘₐₓ)';rr.append(tt);math.append(rr);pp._p.append(math);continue
 if l.startswith('|'):
  group=[l]
  while i<len(lines) and lines[i].startswith('|'):group.append(lines[i]);i+=1
  table(group);continue
 if l.startswith('#'):
  m=re.match(r'(#+) (.*)',l);level=len(m[1]);text=m[2]
  if text.strip().lower()=='references':in_references=True
  if level==1:pp=doc.add_paragraph(text.replace(':','').replace('-',' '),style='Title')
  else:
   numbered=re.match(r'^(\d+(?:\.\d+)*)\.\s+(.*)$',text)
   heading=(numbered.group(1)+'  '+numbered.group(2)) if numbered else text
   heading=re.sub(r'(?<!\d)[/:–—-]|[/:–—-](?!\d)',' ',heading)
   pp=doc.add_paragraph(heading,style=f'Heading {min(level-1,3)}')
  continue
 pp=doc.add_paragraph();inline(pp,l.strip('*') if l.startswith('*') and l.endswith('*') else l)
# Prevent a lone heading but permit natural body pagination.
for pp in doc.paragraphs:pp.paragraph_format.widow_control=True
a.output.parent.mkdir(parents=True,exist_ok=True);doc.save(a.output);print(a.output)
