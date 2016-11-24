classdef OptionStruct < matlab.mixin.Copyable
% Class for option parsing based on structs.
% Construct with OptionStruct(options), where options is either:
%       a list of possible options
%       key value pairs (listing option names and their defaults)
%       a struct
%       an OptionStruct (copies the option struct)
%
% Existing options can be accessed and changed using struct syntax.
% (i.e. OptionStruct.option)
% New options can only be added at construction time.
%
% methods:
%   set('key',value,...): set option values using key-value list
%       (does not create new fields and errors for non-existing
%       options)
%   set(struct): set options using struct (errors for non-existing
%       options)
%
%   setvalid('key',value,...): set valid options from key-value list,
%       returns the unused key-value pairs.
%   
%   setvalid(struct): set valid options from struct, returns the unused
%       key-value pairs as a struct
%
%   isfield('field'): check whether field 'field' exists
%
%   isset('field'): check whether option 'field' is set (i.e. not
%       empty). Can also take a cell array of strings, and than
%       returns a logical vector.
%
%   list(): return options as key-value list
%
%   struct(): return options as struct
%
% properties: options: used to get a list of options
% 
% Version: v1.0-alpha1
% Date: Mon 10 Oct 2016 16:12:35 EDT
% Author: Lucas Jeub
% Email: ljeub@iu.edu
    
    properties (Hidden)
        option_struct=struct([]);
    end
    
    properties (Dependent=true, SetAccess=private)
        options
    end
    
    methods
        
        function obj=OptionStruct(varargin)
            %constructor (pass struct with defaults or list of options)
            if nargin>0
                %unpack nested cells (usually created by passing varargins)
                input=unpack_nested(varargin);
                if length(input)==1
                    if isstruct(input)
                        obj.options=fieldnames(input);
                        obj.set(input);
                    elseif isa(input,'OptionStruct')
                        obj=copy(input);
                    else
                        obj.options=input;
                    end
                else
                    if iscellstr(input)
                        %list of options
                        obj.options=input;
                    elseif ischar(input)
                        %single option
                        ischar(input)
                        obj.options={input};
                    else
                        if ~mod(length(input),2)
                            for i=1:2:length(input)
                                if isvarname(lower(input{i}))
                                    obj.options=lower(input{i});
                                    obj.option_struct.(lower(input{i}))=input{i+1};
                                else
                                    error('%s is not a valid option name',input{i})
                                end
                            end
                        else
                            error('option list has odd length')
                        end
                    end
                end
            end
        end
        
        
        
        function obj=subsasgn(obj, S, value)
            %subscripted assignment
            if isequal(S(1).type,'.')
                if ismember(S(1).subs,properties(obj))||ismember(S(1).subs,methods(obj))
                    obj=builtin('subsasgn',obj,S,value);
                    return;
                end
                
                S(1).subs=lower(S(1).subs);
                if obj.isfield(S(1).subs)
                    obj.option_struct=builtin('subsasgn',obj.option_struct,S,value);
                else
                    error('option %s does not exist',S.subs);
                end
            else
                error('subscripted assignment %s undefined',S.type);
            end
        end
        
        
        function varargout=subsref(obj,S)
            %subscripted reference
            if isequal(S(1).type,'.')
                try
                    if nargout>0
                        varargout{:}=builtin('subsref',obj,S);
                    else
                        builtin('subsref',obj,S);
                    end
                    return;
                catch err
                    if strcmp(err.identifier,'MATLAB:noSuchMethodOrField')
                        S(1).subs=lower(S(1).subs);
                        if obj.isfield(S(1).subs)
                            if nargout>0
                                varargout{:}=builtin('subsref',obj.option_struct,S);
                            else
                                builtin('subsref',obj.option_struct,S);
                            end
                        else
                            error('option %s does not exist',S.subs);
                        end
                    else
                        rethrow(err);
                    end
                end
            else
                error('subscripted reference %s is undefined',S(1).type);
            end
        end
        
        
        function disp(obj)
            %nice display of options
            disp('option struct with fields:')
            disp(obj.option_struct)
        end
        
        %set and get possible option names
        function set.options(obj,fieldnames)
            if ischar(fieldnames)
                fieldnames={fieldnames};
            end
            if ~iscellstr(fieldnames)
                error('cannot parse fieldnames');
            end
            fieldnames=lower(fieldnames);
            for i=1:length(fieldnames)
                obj.option_struct(1).(fieldnames{i})=[];
            end
        end
        function opt=get.options(obj)
            opt=fieldnames(obj.option_struct);
        end
        
        
        function set(obj,varargin)
            % set options (takes a 'struct' or 'key','value' pairs, errors 
            % for bad input)
            input=unpack_nested(varargin);
            if isa(input,'OptionStruct')
                input=input.struct;
            end
            if length(input)==1
                if isstruct(input)
                    fields=fieldnames(input);
                    for i=1:length(fields)
                        if obj.isfield(fields{i})
                            obj.option_struct.(lower(fields{i}))=input.(fields{i});
                        else
                            error('option ''%s'' does not exist',fields{i});
                        end
                    end
                else
                    error('need input struct');
                end
            else
                obj.parse_option_list(input);
            end
        end
        
        function remaining_options=setvalid(obj,varargin)
            % set valid options from input and return the remaining options
            input=unpack_nested(varargin);
            
            if length(input)==1
                if isa(input,'OptionStruct')
                    input=input.struct;
                    return_type='ostruct';
                    fields=fieldnames(input);
                elseif isstruct(input)
                    fields=fieldnames(input);
                    return_type='struct';
                else
                    error('need struct or option struct as input');
                end
            else
                if mod(length(input),2)==0
                    fields=input(1:2:end);
                    values=input(2:2:end);
                    input=cell2struct(values(:),fields,1);
                    return_type='list';
                else
                    error('option list has odd length');
                end
            end
            for i=1:length(fields)
                if obj.isfield(lower(fields{i}))
                    obj.option_struct.(lower(fields{i}))=input.(fields{i});
                    input=rmfield(input,fields{i});
                end
            end
            switch return_type
                case 'ostruct'
                    remaining_options=OptionStruct(input);
                    
                case 'struct'
                    remaining_options=input;
                    
                case 'list'
                    fields=fieldnames(input);
                    values=struct2cell(input);
                    remaining_options=[fields(:)';values(:)'];
                    remaining_options=remaining_options(:)';
            end
        end
        function is_opt=isfield(obj,fieldname)
            % check if fieldname is a valid option
            is_opt=isfield(obj.option_struct,lower(fieldname));
        end
        
        function is_set=isset(obj,fieldname)
            % check if fieldname is set (i.e. not empty)
            fieldname=lower(fieldname);
            if ~iscell(fieldname)
                fieldname={fieldname};
            end
            opt=obj.isfield(fieldname);
            if any(~opt)
                error('option ''%s'' does not exist \n',fieldname{~opt});
            end
            
            is_set=false(size(fieldname));
            for i=1:length(fieldname)
                is_set(i)=~isempty(obj.option_struct.(fieldname{i}));
            end
        end
        
        function options=struct(obj)
            % return options as struct
            options=obj.option_struct;
        end
        
        function options=list(obj)
            % return options as list of key-value pairs
            opts=obj.options;
            vals=struct2cell(obj.option_struct);
            options=[opts(:)';vals(:)'];
            options=options(:)';
        end
        
    end
    
    methods (Access=private)
        function parse_option_list(obj,list)
            if iscell(list)
                if ~mod(length(list),2)
                    for i=1:2:length(list)
                        if obj.isfield(list{i})
                            obj.option_struct.(lower(list{i}))=list{i+1};
                        else
                            optstr=sprintf('%s \n',obj.options{:});
                            error('option ''%s'' does not exist, valid options are %s',list{i},optstr);
                        end
                    end
                else
                    error('option list has odd length')
                end
            else
                error('option list has to be a cell array')
            end
        end
        
        
    end
    
end

function input=unpack_nested(input)
while iscell(input)&&length(input)==1
    input=input{1};
end
end

