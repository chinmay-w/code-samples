// This code snippet is taken from my final project for a Graphics programming class, for which my team chose to implement a "toon" grass shader in Unity HLSL.
// Specifically, it was written to implement shadow casting. While it took me a long time to determine how to best handle the math involved,
// the implementation turned out to be really elegant, for both static and dynamic shadows.

// I have cut out a lot of unrelated code and boilerplate pragra/include directives

// Render pass for our grass shader
Pass
{
    Name "GrsPass"
    Tags {"LightMode" = "UniversalForward"} // to specify Universal Rendering

    HLSLPROGRAM
        float4 fragShader(GeoData g) : SV_Target {

            float4 shadowColor = 1.0f;
            #ifdef _MAIN_LIGHT_SHADOWS // global toggle for shadow rendering in the scene
                VertexPositionInputs input  = (VertexPositionInputs)0;
                input.positionWS = g.worldPos; // Struct that contains the position of point g

                float4 shadowCoord = GetShadowCoord(input); 
                // converts from position into "shadow space" within the Unity library's shadow map
                // The shadow map is made by ray-tracing from all points to the main light source
                // If there is a collision with another object, then that point is in that object's shadow


                float ambient = 0.25f; // Constant used to simulate ambient lighting of a scene
                half attenuation = saturate(MainLightRealtimeShadow(shadowCoord) + ambient); 
                // saturate function wrapper used to clamp between 0 and 1 -- if fully visible then it could have a value of 1.25 for example
                shadowColor = lerp(0.0f, 1.0f, attenuation); // 
            #endif

            // lerp is used to interpolate between base and tip colors for the gradient
            // multplying by the shadow values allows for objects to cast shadows dynamically as they move through a scene
            return shadowColor * lerp(_Base_Color, _Tip_Color, g.uv.y);
        }
    ENDHLSL
}

Pass
{
    name "ShadowCaster"
    {
        "LightMode" = "ShadowCaster" // ShaderLab tag for a pass that renders object depth (from lights into the shadow map)
    }

    HLSLPROGRAM
        float4 fragShaderShadow(GeoData g) : SV_Target {
            // Albedo is a measure of how reflective a surface is, and Alpha is a measure of how opaque a surface is
            // Since both affect shadows cast in the scene, this pass accounts for that when constructing of the shadow map 
            // _BaseMap and _Base_Color are the texture-map and color of the material

            // Shadow maps are much more performant than dynamic shadow simulation, so they are preferred when an object doesn't move
            Alpha(SampleAlbedoAlpha(g.uv, TEXTURE2D_ARGS(_BaseMap, sampler_BaseMap)).a, _Base_Color, 1.0);
            return 0;
        }
    ENDHLSL
}