
from bokeh.models.glyphs import Line
from bokeh.models.selections import IntersectRenderers
from bokeh.palettes import Viridis256
from bokeh.plotting import figure
from bokeh.transform import  factor_cmap
from bokeh.models import LinearColorMapper
import numpy as np


def set_color(df, col, palette = Viridis256):
    """
    Sets the color settings tu use in Cross section plot
    datramapper bases on parameter

    """
    if df[col].apply(isinstance,args = [float]).all():
        mapper = LinearColorMapper(palette=Viridis256, 
                                   low=df[col].min(), 
                                   high=df[col].max())

        color = {'field': col, 'transform': mapper}

    else:

        categories = df[col].drop_duplicates().values

        if len(categories)==1:
            color_list = palette[0:1]

        else:
            color_list = palette[::(len(palette)-1)//(len(categories)-1)]

        color = factor_cmap(col,
                            palette=color_list, 
                            factors=categories)
    return color


def set_str_tooltips(df, *col):
    
    tt_text = [('cell','@ElName')]

    col_unit = {'P':'bar',
                'Pressure':'bar',
                'T':'degC',
                'Temperature':'degC',
                'S_hyd':'frac.',
                'S_Hydrate':'frac.',
                'S_gas':'frac.',
                'S_aqu':'frac.',
                'S_aqueous':'frac.',
                'S_icd':'frac.',
                'S_ice':'frac.',
                'X_inh':'frac.',
                'X_Inhibitor':'frac.',
                'C_CH4inGas':'Kg/m3',
                'C_CH4inAqu':'Kg/m3',
                'Dens_Gas':'Kg/m3',
                'Dens_Aqu':'Kg/m3',
                'Dens_Hydr':'Kg/m3',
                'Visc_Gas':'Kg/m*s',
                'Visc_Aqu':'Kg/m*s',
                'k_rg':'frac.',
                'Krel_Gas':'frac.',
                'k_rw':'frac.',
                'Krel_Aqu':'frac.',
                'k_adj_F':'frac.',
                'perm_abs':'mD',
                'porosity':'frac.',
                'P_CH4':'bar',
                'P_EqHydr':'bar',
                'P_SatWat':'bar',
                'P_cap':'bar',
                'P_aqu':'bar'}
    
    for c in col:
        col_name = c

        if df[c].apply(isinstance,args = [float]).all():

            col_value = '@'+c+' '+col_unit[c]

            tt_text.append((col_name,col_value))

        else:

            col_value = '@'+c


            tt_text.append((col_name,col_value))
    return tt_text
        
def create_figure1(df,col,time_steps,data_source, plane = 'xy'):

    #Define ColumnDataSource
    source_2D=data_source
    
    #Set X, Y ranges
    x=plane[0]
    y=plane[1]

    first_row=df[df[y]==df[y].max()].loc[time_steps[0]]
    last_row= df[df[y]==df[y].min()].loc[time_steps[0]]

    first_col= df[df[x]==df[x].min()].loc[time_steps[0]]
    last_col=df[df[x]==df[x].max()].loc[time_steps[0]]

    Xrange=(first_col[x].min()-first_col['d'+x].max()/2,last_col[x].max()+last_col['d'+x].max()/2)
    Yrange=(last_row[y].min()-last_row['d'+y].max()/2,first_row[y].max()+first_row['d'+y].max()/2)

    #Define
    
    plot_options = dict(toolbar_location=None,
                        y_range=Yrange,
                        x_range=Xrange)    


    figure1 = figure(width=500, x_axis_location='above', plot_height=350,
                     **plot_options)

    label_X = x+' distance [m]'
    label_Y = y+' distance [m]'


    figure1.xaxis.axis_label = label_X

    figure1.yaxis.axis_label = label_Y
    
    figure1.rect(x=x,y=y,width='d'+x,height='d'+y,source=source_2D,
         fill_color=set_color(df,col), line_color=None, line_width=0)  

    return figure1

def create_figure2(df, time_steps, data_source, plane ='xy' ):

    #Define ColumnDataSource
    source_2D=data_source

    #Set X, Y ranges
    x=plane[0]
    y=plane[1]

    first_row=df[df[y]==df[y].max()].loc[time_steps[0]]
    last_row= df[df[y]==df[y].min()].loc[time_steps[0]]

    first_col= df[df[x]==df[x].min()].loc[time_steps[0]]
    last_col=df[df[x]==df[x].max()].loc[time_steps[0]]

    Xrange=(first_col[x].min()-first_col['d'+x].max()/2,last_col[x].max()+last_col['d'+x].max()/2)
    Yrange=(last_row[y].min()-last_row['d'+y].max()/2,first_row[y].max()+first_row['d'+y].max()/2)


    plot_options = dict(toolbar_location=None,
                        y_range=Yrange,
                        x_range=Xrange)    


    figure2 = figure(width=500, x_axis_location='below', plot_height=350,
                    #  tooltips = TOOLTIPS,**plot_options)
                    **plot_options)

    label_X = x+' distance [m]'

    if y == 'z':
        label_Y = 'depth [m]'

    else:
        label_Y = y+' distance [m]'



    figure2.xaxis.axis_label = label_X

    figure2.yaxis.axis_label = label_Y

    #Draw brine
    figure2.rect(x=x,y=y,width='d'+x,height='d'+y,source=source_2D,
        fill_color="#AB957C",fill_alpha={'field':'S_aqu'}, line_color="#AB957C", line_width = 0.5, line_alpha={'field':'S_aqu'}) 

    #Draw gas
    figure2.rect(x=x,y=y,width='d'+x,height='d'+y,source=source_2D,
        fill_color="firebrick",fill_alpha={'field':'S_gas'},line_color="firebrick", line_width = 0.5, line_alpha={'field':'S_gas'}) 



    
    return figure2

def create_figure1_1D(df,col,time_steps,data_source, axis = 'y'):
    #Define ColumnDataSource
    source_1D=data_source

    #Set X, Y ranges
    x=col
    y=axis

    first_row=df[df[y]==df[y].max()].loc[time_steps[0]]
    last_row= df[df[y]==df[y].min()].loc[time_steps[0]]

    if df[col].apply(isinstance,args = [float]).all():
        minv = df[col].min()
        maxv = df[col].max()

        if minv/maxv > 0.97:
            delta = (10**np.floor(np.log10(np.abs(minv))))
            minv -= delta
            maxv += delta

        Xrange= (minv, maxv)

    else:
        Xrange = (0,1)

    Yrange=(last_row[y].min()-last_row['d'+y].max()/2,first_row[y].max()+first_row['d'+y].max()/2)

    #Define

    plot_options = dict(toolbar_location=None,
                        y_range=Yrange, x_range=Xrange)    


    figure1_1D = figure(width=170, x_axis_location='above', y_axis_location='right', plot_height=350,
                    **plot_options)

    label_X = x
    label_Y = axis+' distance [m]'


    figure1_1D.xaxis.axis_label = label_X

    figure1_1D.yaxis.axis_label = label_Y

    figure1_1D.line(x=x, y=y,  line_width=2, color="indigo",source=source_1D)    

    return figure1_1D

def create_figure2_1D(df,time_steps,data_source, axis = 'y'):
    #Define ColumnDataSource
    source_1D=data_source

    #Set X, Y ranges
    y=axis

    first_row=df[df[y]==df[y].max()].loc[time_steps[0]]
    last_row= df[df[y]==df[y].min()].loc[time_steps[0]]


    Yrange=(last_row[y].min()-last_row['d'+y].max()/2,first_row[y].max()+first_row['d'+y].max()/2)

    #Define

    plot_options = dict(toolbar_location=None,
                        y_range=Yrange)    


    figure2_1D = figure(width=170, x_axis_location='below', y_axis_location='right', plot_height=350,
                    **plot_options)

    label_X = 'saturation [frac.]'
    label_Y = axis+' distance [m]'


    figure2_1D.xaxis.axis_label = label_X

    figure2_1D.yaxis.axis_label = label_Y
    #Draw gas

    figure2_1D.line(x='S_aqu', y=y,  line_width=2, color="steelblue",source=source_1D)    
    figure2_1D.line(x='S_gas', y=y,  line_width=2, color="firebrick",source=source_1D)    
    figure2_1D.line(x='S_hyd', y=y,  line_width=2, color="black",source=source_1D)    

    return figure2_1D


def tough_app(doc):
    
    #Retrieve time steps
    tsteps=data.index.get_level_values(0).drop_duplicates()
    tsteps=np.sort(tsteps)
    n_time_steps=len(tsteps)
    
    
    #Selector
    i = 0
    param ='MA12'
    param_select = list(data)[10:]
    
    i_idx_0 = data.I.max()//2
    

     
    #widgets
    select_pmt=Select(value=param, title='Parameter',options=param_select)
    
    slider_time= Slider(start=0, end=n_time_steps-1, value=i, step=1, title="Time step")
    
    slider_1D = Slider(start=data.I.min(), end=data.I.max(), value=i_idx_0, step=1, title="I-index")
    
    i_idx = slider_1D.value

#     slider_X_axis = Slider(start=data.I.min(), end=data.I.max(), value=data.I.max()//2, step=1, title="I-index")
#     slider_Y_axis = Slider(start=data.J.min(), end=data.J.max(), value=data.J.max()//2, step=1, title="J-index")
    
    button_time_next = Button(label="Next")
    button_time_previous = Button(label="Previous")
    button_play = Button(label='► Play', width=80)
    
    
    #Define Data Sources
    data_source_2D=ColumnDataSource(data=data.loc[tsteps[i]])
    data_source_1D=ColumnDataSource(data=data[data.I == slider_1D.value].loc[tsteps[i]])
    
    #Data Source that only contain water column cells
    data_source_wcol = ColumnDataSource(data=data.loc[(data.MA12=='WATER') & (data.index == tsteps[i])])

    
    #Create figures
    top_2D_figure = create_figure1(data,select_pmt.value,tsteps,data_source_2D, plane='xz')
    
    bottom_2D_figure = create_figure2(data,tsteps,data_source_2D, plane='xz')
    
    #Draw water column on 2D figures
    top_2D_figure.rect(x='x', y ='z', width = 'dx', height = 'dz', source = data_source_wcol, fill_color="steelblue", line_width=1, line_color='steelblue')
    bottom_2D_figure.rect(x='x', y ='z', width = 'dx', height = 'dz', source = data_source_wcol, fill_color="steelblue", line_width=1, line_color='steelblue')
    
    top_1D_figure = create_figure1_1D(data,select_pmt.value,tsteps,data_source_1D, axis = 'z')
    top_1D_figure.xaxis.formatter.power_limit_low = -1
    top_1D_figure.xaxis.formatter.precision = 1
    top_1D_figure.xaxis.ticker.desired_num_ticks = 3    
    
    
    bottom_1D_figure = create_figure2_1D(data,tsteps,data_source_1D, axis = 'z')
    
       
    top_2D_figure.add_tools(HoverTool(tooltips = set_str_tooltips(data, select_pmt.value)))
    bottom_2D_figure.add_tools(HoverTool(tooltips = set_str_tooltips(data, 'S_hyd','S_aqu','S_gas')))
    
    Xloc=data[data.I== i_idx].x.max()
    vline = Span(location=Xloc, dimension='height', line_color='red', line_width=1)
    top_2D_figure.renderers.extend([vline])
    bottom_2D_figure.renderers.extend([vline])
    
    #PT-diagram
    if 'P_EqHydr' in param_select:
        
        X_range = (0,25)
        Y_range = (data['P'].min(),data['P'].max())
        
        plot_options = dict(toolbar_location=None,
                            x_range=X_range,
                            y_range=Y_range)  
        
        PT_fig = figure(**plot_options, plot_height=700)
        
        label_X = 'Temperature [degC]'
        label_Y = 'Pressure [MPa]'
        
        
        
        PT_fig.xaxis.axis_label = label_X

        PT_fig.yaxis.axis_label = label_Y
        
        PT_fig.line(x='T', y='P', 
                    color = '#ffa600', 
                    line_width=2, 
                    legend_label = 'Pressure',
                    source = data_source_1D)
        
        PT_fig.line(x='T', y='P_EqHydr', 
                    color = '#bc5090', 
                    line_width=2, 
                    legend_label = 'Stability curve',
                    source = data_source_1D)
        
        
    
    
    VAR=select_pmt.value

    
    def cb_slider_1D(attrname, old, new):
        global i_idx
        i_idx=new
        time_step=slider_time.value
        
        data_source_1D.data=data[data.I == i_idx].loc[tsteps[time_step]].to_dict('list')
        
        vline.location=data[data.I== new].x.max()
    
    def update(attr,old,new):
        #Update colors in 2d Xsections
        top_2D_figure.renderers[0].glyph.update(fill_color=set_color(data,select_pmt.value))
        top_1D_figure.renderers[0].glyph.update(x=select_pmt.value)
        
        #Update data in 1D plots
        top_1D_figure.xaxis.axis_label = select_pmt.value
        top_1D_figure.x_range.reset_start
        top_1D_figure.x_range.reset_end
        
        try:
            
            minv = data[select_pmt.value].min()
            maxv = data[select_pmt.value].max()
           
            if minv == 0 and maxv == 0:
                minv -= 0.5
                maxv += 0.5
                
                
           
            elif minv/maxv > 0.97:
                delta = (10**np.floor(np.log10(np.abs(minv))))/2
                minv -= delta
                maxv += delta

                
               
            
            top_1D_figure.x_range.start = minv
            top_1D_figure.x_range.end = maxv
        
        except:
            top_1D_figure.x_range.start = 0
            top_1D_figure.x_range.end = 1
        

        #Update Hovertool
        
        #Modification due to bug updating value only on tooltips
        top_2D_figure.tools.pop(-1)
        parameter_name = select_pmt.value
        tooltips_str = set_str_tooltips(data, parameter_name)
        top_2D_figure.add_tools(HoverTool(tooltips = tooltips_str))
#         top_2D_figure.hover.tooltips = tooltips_str

        
    def cb_next():
        i = slider_time.value
        i = i + 1
        ev = "button next"
        time_step = i % (n_time_steps)
        slider_time.value = time_step
        
#         print_info(ev,time_step)
        
    def cb_previous():
        i = slider_time.value
        i = i - 1
        ev = "button previous"

        time_step= i % (n_time_steps)
        slider_time.value = time_step
        
#         print_info(ev,time_step)
        
        
    def cb_slider(attrname, old, new):
        global i
        time_step = slider_time.value
        i = slider_time.value
        title.text = '<h3>'+ time_string(tsteps[i]) +'</h3>'
        data_source_2D.data=data.loc[tsteps[time_step]].to_dict('list')
        data_source_wcol.data = data.loc[(data.MA12=='WATER') & (data.index == tsteps[i])].to_dict('list')
        
        data_source_1D.data=data[data.I == slider_1D.value].loc[tsteps[time_step]].to_dict('list')

        

    def animate_update():
        time_step = slider_time.value + 1
        if time_step > n_time_steps-1:
            time_step = 0
        slider_time.value = time_step
    
    
    callback_id = None
    
    def animate():
        global callback_id
        if button_play.label == '► Play':
            button_play.label = '❚❚ Pause'
            callback_id = doc.add_periodic_callback(animate_update, 200)
        else:
            button_play.label = '► Play'
            doc.remove_periodic_callback(callback_id)





    select_pmt.on_change('value',update)
    
    slider_time.on_change('value',cb_slider)
    slider_1D.on_change('value',cb_slider_1D)
    button_time_next.on_click(cb_next)
    button_time_previous.on_click(cb_previous)
    button_play.on_click(animate)

    
    buttons_row = row(button_play,button_time_previous, slider_time, button_time_next, slider_1D,width=800)
    top_row = row(top_2D_figure,top_1D_figure )
    bottom_row = row(bottom_2D_figure, bottom_1D_figure)
    title = Div(text='<h3>'+ time_string(tsteps[i]) +'</h3>')
    
    if 'P_EqHydr' in param_select:
        display_layout = layout([title],
                                [select_pmt],
                                [buttons_row],
                                [[[top_row],
                                [bottom_row]],
                                PT_fig])
    
    else:    
        display_layout = layout([title],
                                [select_pmt],
                                [buttons_row],
                                [top_2D_figure,top_1D_figure],
                                [bottom_2D_figure,bottom_1D_figure])

        
    doc.add_root(display_layout)