if JSON and JSONSchemaValidator then
    --print("JSON setup BEGIN")
    --print("* Reading schema...")
    local schema=JSON(); 
    local my_path=debug.getinfo(1).source:match("@?(.*/)")
    JSON_load_from_file(schema, my_path.."json/schema/problem-schema.json");
    -- print(JSON_dump(schema))
    
    --print("* Creating validator...")
    local validator=JSONSchemaValidator();
    validator:init(schema)

    local json = 
    {
        validate = function (luadesc)
            -- print("* Reading problem...")
            local document=JSON(); 
            JSON_parse(document, util.json.encode(luadesc))  
            -- print("* Validating problem...")
            return validator:validate_document(document)
        end    
    }

    --print("JSON setup END")
    return json
else
    print("JSON validation is NOT AVAILABLE.")
    return nil
end