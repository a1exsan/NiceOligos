FROM python:3.9-slim

EXPOSE 8080

WORKDIR /app

#RUN apt-get update && apt-get install -y libxrender1 libxext6 libsm6 && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y libxrender1 libxext6 libsm6 libexpat1 && rm -rf /var/lib/apt/lists/*

ADD requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

ADD external_mods.csv .
ADD chemicals_page.py .
ADD input_order_page.py .
ADD invoce_chart.py .
ADD invoce_page.py .
ADD main.py .
ADD OligoMap_utils.py .
ADD oligosynth_page.py .
ADD raw_material_page.py .
ADD raw_mat_widget.py .
ADD Reactor.py .
ADD xwell_plate_unit.py .
ADD lcms_chrom_data.py .
ADD stock_data_page.py .
ADD lcms_dialog_model.py .
ADD lcms_zip.py .
ADD mzdatapy.py .
ADD server.py .
ADD molseq_lang.py .
ADD synthesis_method.py .
ADD asm2000_method.py .
ADD asm_templates.py .
ADD mol3D.py .

RUN mkdir -p /app/static_images
RUN mkdir -p /app/images
RUN mkdir -p /app/templates

ADD static_images/chart_bkg_1.png ./static_images
ADD static_images/chart_bkg_2.png ./static_images
ADD static_images/chart_bkg_3.png ./static_images
ADD static_images/favicon.ico ./static_images
ADD static_images/icon.png ./static_images
ADD static_images/infopanel_bkg_1.png ./static_images
ADD static_images/infopanel_bkg_2.png ./static_images
ADD static_images/number_oligos_plot_1.png ./static_images
ADD static_images/plate_bkg_1.png ./static_images
ADD static_images/plate_bkg_2.png ./static_images
ADD static_images/widget_bkg_1.png ./static_images
ADD static_images/widget_bkg_2.png ./static_images
ADD static_images/widget_bkg_3.png ./static_images
ADD static_images/widget_bkg_4.png ./static_images

ADD images/background_1.jpeg ./images

ADD templates/pass_tmpl.docx ./templates
ADD templates/pass_tmpl_2.docx ./templates
ADD templates/combined.docx ./templates
ADD templates/passport_doc.docx ./templates
ADD templates/passport_rowdata.docx ./templates
ADD templates/autosampler_75ul_method1.pr2 ./templates
ADD templates/autosampler_75ul_template.pr2 ./templates

ADD templates/lcms_1D_plot.png ./templates
ADD templates/lcms_2D_plot.png ./templates
ADD templates/lcms_report_tmpl_1.docx ./templates
ADD templates/lcms_report_tmpl_2.docx ./templates
ADD templates/lcms_report_tab.docx ./templates
ADD templates/combined_lcms_report.docx ./templates
ADD templates/lcms_report_doc.docx ./templates

# VOLUME "/app/data/temp"

ENTRYPOINT ["python"]
CMD ["main.py"]
