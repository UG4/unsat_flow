local json = require("json")
local template_import = require("config/template")

function problem_check(problem)
    if problem == nil then
        print("problem not defined or not named 'problem'")
        return false
    end

    function obj_in_table(object, table)
        for i, v in pairs(table) do
            if v == object then
                return true
            end
        end
        return false
    end

    function type_check(problem_type, template)
        if type(template) == "string" then
            if type(problem_type) ~= template then
                return false
            end
        elseif type(template) == "table" then
            if  obj_in_table(type(problem_type), template) or
                obj_in_table(problem_type, template) then
                return true
            else
                return false
            end
        end

        return true
    end

    function table_content_type(ctype, table)
        for i, v in ipairs(table) do
            if type(v) ~= ctype then
                return false
            end
        end
        return true
    end

    function value_check(value, template_list)
        for k, v in pairs(template_list) do
            -- there are a few keywords that can be parsed into the lists in the 
            -- template file. current possibilities are:
            -- __min: a minimal value for a number
            -- __max: respectively a max value for a number
            -- __type: a table containing possible types for a variable
            -- __content_type: checks if the content of a table is of this type
            if k == "__min" then
                if value < v then
                    return false
                end
            elseif k == "__max" then
                if value > v then
                    return false
                end
            elseif k == "__type" then
                if not type_check(value, v) then
                    return false
                end
            elseif k == "__content_type" then
                if not table_content_type(v, value) then
                    return false
                end
            end
        end
        return true
    end

    function table_check(table, template)
        if table ~= nil then
            if template["__repeatable"] == true then
                for k, v in ipairs(table) do
                    if not table_check(v, template["__values"]) then
                        return false
                    end
                end
            else
                if template["__type"] ~= nil then
                    if template["__deftypes"] == true then
                        for i, v in ipairs(template["__type"]) do
                            if table.type == v then
                                if not table_check(table, template["__"..v]) then
                                    return false
                                end
                                return true
                            end
                        end
                        print("Could not find a type")
                        return false
                    end

                    if not value_check(table, template) then
                        return false
                    end
                else
                    for key, value in pairs(template) do
                        if not table_check(table[key], value) then
                            print("failed at key ", key)
                            print("value should be ", value)
                            print("actual value is ", table[key])
                            return false
                        end
                    end
                end
            end
        else
            return false
        end
        return true
    end

    for key, value in pairs(template) do
        if not table_check(problem[key], value) then
            return false
        end
    end
    print("problem configuration correct")

    return true
end