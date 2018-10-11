classdef SystemsGUI < handle
    
    properties
        Figure
        
        nameBox
        variablesBox
        parametersBox
        timeBox
        derivativesButtonGroup
        equationsBox
        myWaitbar
        progress = 0        
    end
    
    methods

        % constructor
        function app = SystemsGUI

            pitch = 30;

            screensize = get(0,'ScreenSize');
            screenwidth = screensize(3);
            screenheight = screensize(4);
            windowwidth = 0.5 * screenwidth;
            windowheight = 0.75 * screenheight;
            
            title = 'MatcontL system editor';
            
            
            app.Figure = figure('MenuBar','none',... % Main figure
                'NumberTitle','off',...
                'Name',title,...
                'CloseRequestFcn',@app.closeApp,...
                'OuterPosition',...
                    [0.07*screenwidth,0.18*screenheight,...
                     windowwidth,windowheight],...
                'Resize','off');
            
            innerWindowSize = get(app.Figure,'InnerPosition');
            innerWidth = innerWindowSize(3);
            innerHeight = innerWindowSize(4);
            
            labelwidth = 190;
            app.createLabel(...
                [30 innerHeight - 1.5 *  pitch labelwidth 14],'System name');
            app.createLabel(...
                [30 innerHeight - 2.5 * pitch labelwidth 14],'Variables');
            app.createLabel(...
                [30 innerHeight - 3.5 * pitch labelwidth 14],'Parameters');
            app.createLabel(...
                [30 innerHeight - 4.5 * pitch labelwidth 14],'Time variable');
            app.createLabel(...
                [30 innerHeight - 5.75 * pitch labelwidth 14],...
                'Maximum order of symbolic derivatives');
            app.createLabel(...
                [30 innerHeight - 6.75* pitch labelwidth 14],'System equations');
            
            app.nameBox = app.createTextbox(...
                [30+labelwidth+30, innerHeight - 1.5 * pitch,...
                innerWidth - (30+labelwidth+30) - 30, 20]);
            app.variablesBox = app.createTextbox(...
                [30+labelwidth+30, innerHeight - 2.5 * pitch,...
                innerWidth - (30+labelwidth+30) - 30, 20]);
            app.parametersBox = app.createTextbox(...
                [30+labelwidth+30, innerHeight - 3.5 * pitch,...
                innerWidth - (30+labelwidth+30) - 30, 20]);
            app.timeBox = app.createTextbox(...
                [30+labelwidth+30, innerHeight - 4.5 * pitch,...
                innerWidth - (30+labelwidth+30) - 30, 20]);
            
            app.timeBox.String = "t";
            
            app.derivativesButtonGroup = uibuttongroup(app.Figure,...
                'Visible','off',...
                'Units','pixels',...
                'Position', [30+labelwidth+30, innerHeight - 6 * pitch,...
                innerWidth - (30+labelwidth+30) - 30, 30]...
               );
            
            % note: coordinates of the radiobuttons
            % are relative to the uibuttongroup
            app.createRadiobutton([20     ,0,50,30],'1');
            app.createRadiobutton([20+1*50,0,50,30],'2');
            app.createRadiobutton([20+2*50,0,50,30],'3');
            app.createRadiobutton([20+3*50,0,50,30],'4');
            app.createRadiobutton([20+4*50,0,50,30],'5');
           
            app.equationsBox = uicontrol(app.Figure,...
                'Style','edit','Position',...
                [30 60 innerWidth - 2*30 innerHeight - 9 * pitch],...
                'Max',2,...
                'HorizontalAlignment','left');
            
            app.derivativesButtonGroup.Visible = 'on';
            
            uicontrol(app.Figure,'Style','pushbutton','String','OK',...
                'Position',[ innerWidth-245 15 100 30],...
                'Callback',@app.OKCallback);
            uicontrol(app.Figure,'Style','pushbutton','String','Cancel',...
                'Position',[ innerWidth-130 15 100 30],...
                'Callback',@app.closeApp);
        end
 
        function label=createLabel(app,position,text)
            label = uicontrol(app.Figure,...
                'Style','text','Position',position,...
                'HorizontalAlignment','left',...
                'String',text);
        end
        
        function textbox=createTextbox(app,position)
            textbox = uicontrol(app.Figure,...
                'Style','edit','Position',position,...
                'HorizontalAlignment','left');
        end
        
        function textbox=createRadiobutton(app,position,text)
            textbox = uicontrol(app.derivativesButtonGroup,...
                'Style','radiobutton','Position',position,...
                'String',text);
        end

        
    function OKCallback(app,~,~)

      name = app.nameBox.String;
      variables =  app. variablesBox.String;
      parameters = app.parametersBox.String;
      time =             app.timeBox.String;
      maxOrder = str2num( ...
        app.derivativesButtonGroup.SelectedObject.String); %#ok<ST2NM>
      equations =  app.equationsBox.String;
      try
        app.progress = 0;
        app.myWaitbar = waitbar(app.progress, 'Please wait...');
        s = System_of_ODEs(name,variables,parameters,time, ...
          maxOrder,equations,app);
        s.generate_file
        close(app.myWaitbar)
        msgbox(['System file ' char(s.name) '.m was created']);
      catch ex
        if   strcmp(ex.identifier,'System_of_ODEs:BadInput')...
          || strcmp(ex.identifier,'System_of_ODEs:EvalError')...
          || strcmp(ex.identifier,'System_of_ODEs:UndefinedFunctionInEval')
                    errordlg(ex.message, 'MatcontL Error');
        else
          ex.rethrow
        end
      end
    end
    
    function updateStatus(app,str)
      app.progress = app.progress + 0.2;
      waitbar(app.progress, app.myWaitbar, str);
    end
        
    function closeApp(app,~,~)
      delete(app.Figure)
    end

  end
end                